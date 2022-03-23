function res=runNWT(nwt)

%Save specification
res.nwtSpec=nwt;

%Perform time integration
tic
switch nwt.solver.type
    case 'hos'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dx=Lx/Nx;
                x=(dx*(0:(Nx-1)));
                
                if nwt.sim.dt ~= 0
                    returnTimes = 0:nwt.sim.dt:nwt.sim.tMax;
                else
                    returnTimes = [0,nwt.sim.tMax];
                end
                
                %Initialize surface elevation and potential with linear solution
                eta0=zeros(Nx,1);
                phiS0=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            if isfield(nwt.init.waveComp{iComp},'fHigh')
                                if isfield(nwt.init.waveComp{iComp},'filter')
                                    switch nwt.init.waveComp{iComp}.filter.type
                                        case 'power'
                                            kHigh=level1.computeWaveNumber(2*pi*nwt.init.waveComp{iComp}.fHigh,h,0,0);
                                            HLPJ=1./(1+(k/kHigh).^nwt.init.waveComp{iComp}.filter.n).';
                                            psiW=psiW.*HLPJ;
                                            psiW(psiW<1e-8)=0;
                                            psiW=psiW*(Hs/(4*sqrt(sum(psiW)*dk)))^2;
                                        otherwise
                                            error('a:a','Unknown filter type %s\n',nwt.init.waveComp{iComp}.filter.type)
                                    end
                                else
                                    error('a:a','Missing filter definition for wave component #%d\n',iComp)
                                end
                            end
                            [etaTemp,phiSTemp]=initializeEtaPhiS_Irregular(Nx,h,k.',psiW,seed);
                            eta0=eta0+etaTemp;
                            phiS0=phiS0+phiSTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [etaTemp,phiSTemp]=initializeEtaPhiS_ComplexAmplitudes(k.',compAmp,h);
                            eta0=eta0+etaTemp;
                            phiS0=phiS0+phiSTemp;
                        case 'Gaussian'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            sigmak=nwt.init.waveComp{iComp}.sigmak;
                            k0=nwt.init.waveComp{iComp}.k0;
                            seed=nwt.init.waveComp{iComp}.seed;
                            psiW=computeWaveSpectralDensity_k_Gaussian(k.',Hs,k0,sigmak);
                            [etaTemp,phiSTemp]=initializeEtaPhiS_Irregular(Nx,h,k.',psiW,seed);
                            eta0=eta0+etaTemp;
                            phiS0=phiS0+phiSTemp;
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[eta0;phiS0];
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=dk*nwt.solver.LPfilter.kCutMode;
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    case 'cut'
                        HLP=ones(size(k));
                        HLP(nwt.solver.LPfilter.kCutMode+1:end)=0;
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        nRamp=nwt.solver.ramp.n;
                        [t,Y]=ode45(@(t,y)detaPhiS_dt_ExpRamp(t,y,k,M,h,HLP,TRamp,nRamp),returnTimes,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.x=x;
        res.eta=Y(:,1:Nx);
        res.phiS=Y(:,Nx+1:2*Nx);
        res.simTime=toc;
        
    case 'melhos'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1)));
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                deltax0=zeros(Nx,1);
                deltaz0=zeros(Nx,1);
                phiLS0=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            if isfield(nwt.init.waveComp{iComp},'fHigh')
                                if isfield(nwt.init.waveComp{iComp},'filter')
                                    switch nwt.init.waveComp{iComp}.filter.type
                                        case 'power'
                                            kHigh=level1.computeWaveNumber(2*pi*nwt.init.waveComp{iComp}.fHigh,h,0,0);
                                            HLPJ=1./(1+(k/kHigh).^nwt.init.waveComp{iComp}.filter.n).';
                                            psiW=psiW.*HLPJ;
                                            psiW(psiW<1e-8)=0;
                                            psiW=psiW*(Hs/(4*sqrt(sum(psiW)*dk)))^2;
                                        otherwise
                                            error('a:a','Unknown filter type %s\n',nwt.init.waveComp{iComp}.filter.type)
                                    end
                                else
                                    error('a:a','Missing filter definition for wave component #%d\n',iComp)
                                end
                            end
                            [deltaxTemp,deltazTemp,phiLSTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_ComplexAmplitudes(k.',compAmp,h);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Gaussian'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            sigmak=nwt.init.waveComp{iComp}.sigmak;
                            k0=nwt.init.waveComp{iComp}.k0;
                            seed=nwt.init.waveComp{iComp}.seed;
                            psiW=computeWaveSpectralDensity_k_Gaussian(k.',Hs,k0,sigmak);
                            [deltaxTemp,deltazTemp,phiLSTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Regular_Order5'
                            k0=dk*nwt.init.waveComp{iComp}.nk0;
                            a=nwt.init.waveComp{iComp}.a;
                            [deltax0,deltaz0,phiLS0]=initializeXZPhiLS_Order5(k0,a,h,alpha);
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[deltax0;deltaz0;phiLS0];
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=dk*nwt.solver.LPfilter.kCutMode;
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    case 'cut'
                        HLP=ones(size(k));
                        HLP(nwt.solver.LPfilter.kCutMode+1:end)=0;
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        nRamp=nwt.solver.ramp.n;
                        if isfield(nwt.solver,'noiseRemoval')
                            tStep=nwt.solver.noiseRemoval.resetTime;
                            tStart=0;
                            tEnd=min([tStep,tMax]);
                            isOneMore=false;
                            t=[];
                            Y=[];
                            while (tEnd<=tMax)||isOneMore
                                [ti,Yi]=ode45(@(t,y)dXZPhiLS_dt_ExpRamp(t,y,k,M,h,HLP,TRamp,nRamp),tStart:dt:tEnd,y0,options);
                                fprintf('%8.3f - %8.3f\n',ti(1),ti(end));
                                tStart=tEnd;
                                tEnd=tEnd+tStep;
                                xi=Yi(end,1:Nx);
                                zi=Yi(end,Nx+1:2*Nx);
                                phiLSi=Yi(end,2*Nx+1:end);
                                Axi=getModeAmplitudes(xi);
                                Azi=getModeAmplitudes(zi);
                                AphiLSi=getModeAmplitudes(phiLSi);
                                figure;plot(log10(1e-10+abs(Axi)));grid on;title(['t=',num2str(ti(end),'%.2f')])
                                Axi(abs(Axi)<1e-6)=0;
                                Azi(abs(Azi)<1e-6)=0;
                                AphiLSi(abs(AphiLSi)<1e-6)=0;
                                y0=[getFunctionFromModeAmplitudes(Axi),getFunctionFromModeAmplitudes(Azi),getFunctionFromModeAmplitudes(AphiLSi)].';
                                if tEnd<=tMax&&~isOneMore
                                    isOneMore=true;
                                else
                                    isOneMore=false;
                                end
                                if ~isOneMore
                                    t=[t;ti]; %#ok<AGROW>
                                    Y=[Y;Yi]; %#ok<AGROW>
                                else
                                    t=[t;ti(1:end-1)]; %#ok<AGROW>
                                    Y=[Y;Yi(1:end-1,:)]; %#ok<AGROW>
                                end
                            end
                        else
                            [t,Y]=ode45(@(t,y)dXZPhiLS_dt_ExpRamp(t,y,k,M,h,HLP,TRamp,nRamp),0:dt:tMax,y0,options);
                        end
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=Y(:,1:Nx);
        res.z=Y(:,Nx+1:2*Nx);
        res.phiLS=Y(:,2*Nx+1:end);
        %         res.etaRecons=zeros(size(res.z));
        res.eta=zeros(size(res.z));
        for it=1:length(t)
            deltaXPrime=mean(res.x(it,:));
            deltaX=res.x(it,:)-deltaXPrime;
            etaPrime=sum(computeEta(deltaX,res.z(it,:),k,M),1);
            res.z(it,:)=res.z(it,:)-mean(etaPrime);
            [alphaSorted,indSort]=sort(mod(alpha+deltaXPrime,nwt.sim.Lx));
            %             etaSorted=etaPrime(indSort);
            xSorted=res.x(it,indSort)-deltaXPrime;
            zSorted=res.z(it,indSort);
            res.eta(it,:)=interp1([alphaSorted(end-Nx/4:end)+xSorted(end-Nx/4:end)-nwt.sim.Lx,alphaSorted+xSorted,alphaSorted(1:Nx/4)+xSorted(1:Nx/4)+nwt.sim.Lx],...
                [zSorted(end-Nx/4:end),zSorted,zSorted(1:Nx/4)],alpha);
            %             res.etaRecons(it,:)=interp1([alphaSorted(end)-nwt.sim.Lx,alphaSorted,alphaSorted(1)+nwt.sim.Lx],...
            %                 [etaSorted(end),etaSorted,etaSorted(1)],alpha)-mean(etaPrime);
        end
        res.simTime=toc;
        
    case 'melhos_PressCorr'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1)));
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                deltax0=zeros(Nx,1);
                deltaz0=zeros(Nx,1);
                phiLS0=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            psiW(k>2*pi/nwt.init.lambdaCut)=0;
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_ComplexAmplitudes(k.',compAmp,h);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Gaussian'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            sigmak=nwt.init.waveComp{iComp}.sigmak;
                            k0=nwt.init.waveComp{iComp}.k0;
                            seed=nwt.init.waveComp{iComp}.seed;
                            psiW=computeWaveSpectralDensity_k_Gaussian(k.',Hs,k0,sigmak);
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[deltax0;deltaz0;phiLS0];
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=k(nwt.solver.LPfilter.kCutMode);
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Precompute abcd-coefficients for convective acceleration
                abcd=getabcd(M);
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        nRamp=nwt.solver.ramp.n;
                        [t,Y]=ode45(@(t,y)dXZPhiLS_dt_ExpRampPressCorr(t,y,k,M,h,HLP,TRamp,nRamp,abcd),0:dt:tMax,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=Y(:,1:Nx);
        res.z=Y(:,Nx+1:2*Nx);
        res.phiLS=Y(:,2*Nx+1:end);
        res.eta=zeros(size(res.z));
        res.ax=zeros(size(res.z));
        res.az=zeros(size(res.z));
        res.zCorr=zeros(size(res.z));
        for it=1:length(t)
            deltaXPrime=mean(res.x(it,:));
            deltaX=res.x(it,:)-deltaXPrime;
            AphiS=getModeAmplitudes(res.phiLS(it,:));
            res.zCorr(it,:)=computeZPressCorr(deltaX,res.z(it,:),AphiS,k,M,h,abcd,1);
            try
                [res.ax(it,:),res.az(it,:)]=computeParticleAcc(deltaX,res.z(it,:),AphiS,k,M,h,abcd,1);
            catch
                keyboard
            end
            res.z(it,:)=res.z(it,:)+res.zCorr(it,:);
            [alphaSorted,indSort]=sort(mod(alpha+deltaXPrime,nwt.sim.Lx));
            xSorted=res.x(it,indSort)-deltaXPrime;
            zSorted=res.z(it,indSort);
            res.eta(it,:)=interp1([alphaSorted(end-Nx/4:end)+xSorted(end-Nx/4:end)-nwt.sim.Lx,alphaSorted+xSorted,alphaSorted(1:Nx/4)+xSorted(1:Nx/4)+nwt.sim.Lx],...
                [zSorted(end-Nx/4:end),zSorted,zSorted(1:Nx/4)],alpha);
        end
        res.simTime=toc;
        
    case 'melhos_kInt'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1)));
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                deltax0=zeros(Nx,1);
                deltaz0=zeros(Nx,1);
                phiLS0=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            psiW(k>2*pi/nwt.init.lambdaCut)=0;
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [deltaxTemp,deltazTemp,phiLSTemp,~]=initializeXZPhiLS_ComplexAmplitudes(k.',compAmp,h);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[getModeAmplitudes(deltax0.').';getModeAmplitudes(deltaz0.').';getModeAmplitudes(phiLS0.').'];
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=k(nwt.solver.LPfilter.kCutMode);
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        nRamp=nwt.solver.ramp.n;
                        [t,Y]=ode45(@(t,y)dXZPhiLS_dt_ExpRamp_kInt(t,y,k,M,h,HLP,TRamp,nRamp),0:dt:tMax,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=getFunctionFromModeAmplitudes(Y(:,1:Nx/2));
        res.z=getFunctionFromModeAmplitudes(Y(:,Nx/2+1:2*Nx/2));
        res.phiLS=getFunctionFromModeAmplitudes(Y(:,2*Nx/2+1:end));
        res.eta=zeros(size(res.z));
        for it=1:length(t)
            deltaX=res.x(it,:)-mean(res.x(it,:));
            res.eta(it,:)=sum(computeEta(deltaX,res.z(it,:),k,M),1);
            res.z(it,:)=res.z(it,:)-mean(res.eta(it,:));
        end
        res.simTime=toc;
        
    case 'melhos_acc'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1)));
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                deltax0=zeros(Nx,1);
                deltaz0=zeros(Nx,1);
                phiLS0=zeros(Nx,1);
                dphiLS0dt=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            psiW(k>2*pi/nwt.init.lambdaCut)=0;
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_ComplexAmplitudes(k.',compAmp,h);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        case 'Gaussian'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            sigmak=nwt.init.waveComp{iComp}.sigmak;
                            k0=nwt.init.waveComp{iComp}.k0;
                            seed=nwt.init.waveComp{iComp}.seed;
                            psiW=computeWaveSpectralDensity_k_Gaussian(k.',Hs,k0,sigmak);
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[deltax0;deltaz0;phiLS0;dphiLS0dt;0];
                
                %Precompute abcd-coefficients for convective acceleration term
                abcd=getabcd(M);
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=k(nwt.solver.LPfilter.kCutMode);
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        tau=nwt.solver.relaxTimeMWL;
                        nRamp=nwt.solver.ramp.n;
                        [t,Y]=ode45(@(t,y)dXZPhiLS_dt_AccExpRamp(t,y,k,M,h,HLP,TRamp,nRamp,tau,abcd),0:dt:tMax,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=Y(:,1:Nx);
        res.z=Y(:,Nx+1:2*Nx);
        res.phiLS=Y(:,2*Nx+1:3*Nx);
        res.dphiLS_dt=Y(:,3*Nx+1:4*Nx);
        res.zm=Y(:,end);
        res.z=res.z-repmat(res.zm,1,size(res.z,2));
        res.eta=zeros(size(res.z));
        for it=1:length(t)
            deltaXPrime=mean(res.x(it,:));
            [alphaSorted,indSort]=sort(mod(alpha+deltaXPrime,nwt.sim.Lx));
            xSorted=res.x(it,indSort)-deltaXPrime;
            zSorted=res.z(it,indSort);
            res.eta(it,:)=interp1([alphaSorted(end-Nx/4:end)+xSorted(end-Nx/4:end)-nwt.sim.Lx,alphaSorted+xSorted,alphaSorted(1:Nx/4)+xSorted(1:Nx/4)+nwt.sim.Lx],...
                [zSorted(end-Nx/4:end),zSorted,zSorted(1:Nx/4)],alpha);
        end
        res.simTime=toc;
        
    case 'melhos_accPCorr'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1)));
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                deltax0=zeros(Nx,1);
                deltaz0=zeros(Nx,1);
                phiLS0=zeros(Nx,1);
                dphiLS0dt=zeros(Nx,1);
                for iComp=1:length(nwt.init.waveComp)
                    switch nwt.init.waveComp{iComp}.type
                        case 'JONSWAP'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            Tp=nwt.init.waveComp{iComp}.Tp;
                            seed=nwt.init.waveComp{iComp}.seed;
                            gam=nwt.init.waveComp{iComp}.gamma;
                            psiW=computeWaveSpectralDensity_k_JONSWAP(k.',Hs,Tp,gam,h);
                            psiW(k>2*pi/nwt.init.lambdaCut)=0;
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        case 'Regular_ModeNumber'
                            compAmp=zeros(Nx/2,1);
                            compAmp(nwt.init.waveComp{iComp}.nk0+1)=nwt.init.waveComp{iComp}.a*exp(1i*nwt.init.waveComp{1}.phaseRad);
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_ComplexAmplitudes(k.',compAmp,h);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        case 'Gaussian'
                            Hs=nwt.init.waveComp{iComp}.Hs;
                            sigmak=nwt.init.waveComp{iComp}.sigmak;
                            k0=nwt.init.waveComp{iComp}.k0;
                            seed=nwt.init.waveComp{iComp}.seed;
                            psiW=computeWaveSpectralDensity_k_Gaussian(k.',Hs,k0,sigmak);
                            [deltaxTemp,deltazTemp,phiLSTemp,dphiLSdtTemp]=initializeXZPhiLS_Irregular(Nx,h,k.',psiW,seed);
                            deltax0=deltax0+deltaxTemp;
                            deltaz0=deltaz0+deltazTemp;
                            phiLS0=phiLS0+phiLSTemp;
                            dphiLS0dt=dphiLS0dt+dphiLSdtTemp;
                        otherwise
                            error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{iComp}.type)
                    end
                end
                y0=[deltax0;deltaz0;phiLS0;dphiLS0dt;0];
                
                %Precompute abcd-coefficients for convective acceleration term
                abcd=getabcd(M);
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=k(nwt.solver.LPfilter.kCutMode);
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        tauMWL=nwt.solver.relaxTimeMWL;
                        tauPress=nwt.solver.relaxTimePress;
                        nRamp=nwt.solver.ramp.n;
                        [t,Y]=ode45(@(t,y)dXZPhiLS_dt_AccPCorrExpRamp(t,y,k,M,h,HLP,TRamp,nRamp,tauMWL,tauPress,abcd),0:dt:tMax,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=Y(:,1:Nx);
        res.z=Y(:,Nx+1:2*Nx);
        res.phiLS=Y(:,2*Nx+1:3*Nx);
        res.dphiLS_dt=Y(:,3*Nx+1:4*Nx);
        res.zm=Y(:,end);
        res.z=res.z-repmat(res.zm,1,size(res.z,2));
        res.eta=zeros(size(res.z));
        for it=1:length(t)
            deltaXPrime=mean(res.x(it,:));
            [alphaSorted,indSort]=sort(mod(alpha+deltaXPrime,nwt.sim.Lx));
            xSorted=res.x(it,indSort)-deltaXPrime;
            zSorted=res.z(it,indSort);
            res.eta(it,:)=interp1([alphaSorted(end-Nx/4:end)+xSorted(end-Nx/4:end)-nwt.sim.Lx,alphaSorted+xSorted,alphaSorted(1:Nx/4)+xSorted(1:Nx/4)+nwt.sim.Lx],...
                [zSorted(end-Nx/4:end),zSorted,zSorted(1:Nx/4)],alpha);
        end
        res.simTime=toc;
        
    case 'melhos_Order5'
        M=nwt.solver.M;
        
        switch nwt.sim.type
            case 'periodicDomain2D'
                
                %Domain parameters
                Nx=nwt.sim.Nx;
                Lx=nwt.sim.Lx;
                dk=2*pi/Lx;
                k=dk*(0:(Nx/2-1));
                h=nwt.sim.depth;
                dalpha=Lx/Nx;
                alpha=(dalpha*(0:(Nx-1))).';
                
                %Time integration parameters
                dt=nwt.sim.dt;
                tMax=nwt.sim.tMax;
                
                %Initialize particle displacement and potential with linear solution
                switch nwt.init.waveComp{1}.type
                    case 'Regular_Order5'
                        k0=dk*nwt.init.waveComp{1}.nk0;
                        a=nwt.init.waveComp{1}.a;
                        [deltax0,deltaz0,phiLS0]=initializeXZPhiLS_Order5(k0,a,h,alpha);
                    otherwise
                        error('a:a','Unknown wave component type %s\n',nwt.init.waveComp{1}.type)
                end
                y0=[deltax0;deltaz0;phiLS0];
                
                %Define low-pass filter to remove HF numerical noise
                switch nwt.solver.LPfilter.type
                    case 'power'
                        kCut=dk*nwt.solver.LPfilter.kCutMode;
                        HLP=1./(1+(k/kCut).^nwt.solver.LPfilter.n);
                    otherwise
                        error('a:a','Unknown Low-pass filter type %s\n',nwt.solver.LPfilter.type)
                end
                
                %Perform time-integration
                rTol=nwt.solver.rTol;
                options=odeset('RelTol',rTol);
                switch nwt.solver.ramp.type
                    case 'exp'
                        TRamp=nwt.solver.ramp.Ta;
                        nRamp=nwt.solver.ramp.n;
                            [t,Y]=ode45(@(t,y)dXZPhiLS_dt_ExpRampO5(t,y,k,M,h,HLP,TRamp,nRamp),0:dt:tMax,y0,options);
                    otherwise
                        error('a:a','Unknown ramp type %s\n',nwt.hos.ramp.type)
                end
            otherwise
                error('a:a','Unknown hos simulation type %s\n',nwt.hos.simType)
                
        end
        
        %Extract MELHOS-simulation results + postprocessing of Eulerian surface elevation
        res.t=t;
        res.alpha=alpha;
        res.x=Y(:,1:Nx);
        res.z=Y(:,Nx+1:2*Nx);
        res.phiLS=Y(:,2*Nx+1:end);
        res.eta=zeros(size(res.z));
        for it=1:length(t)
            deltaXPrime=mean(res.x(it,:));
            deltaX=res.x(it,:)-deltaXPrime;
            etaPrime=sum(computeEta(deltaX,res.z(it,:),k,M),1);
            res.z(it,:)=res.z(it,:)-mean(etaPrime);
            [alphaSorted,indSort]=sort(mod(alpha+deltaXPrime,nwt.sim.Lx));
            xSorted=res.x(it,indSort)-deltaXPrime;
            zSorted=res.z(it,indSort);
            res.eta(it,:)=interp1([alphaSorted(end-Nx/4:end)+xSorted(end-Nx/4:end)-nwt.sim.Lx,alphaSorted+xSorted,alphaSorted(1:Nx/4)+xSorted(1:Nx/4)+nwt.sim.Lx],...
                [zSorted(end-Nx/4:end),zSorted,zSorted(1:Nx/4)],alpha);
        end
        res.simTime=toc;
        
    otherwise
        error('a:a','Unknown solver type %s\n',nwt.solver);
end

