clear all

csvFile = 'C:\Users\andreasak\SINTEF\OSC - Current system CFD simulations folder - Dokumenter\BasinDesignSep2021\CFD_results\4Csaba\Case000_UNS\XYZ_Internal_Table_1000.csv';
tab = readtable(csvFile);

assert( numel(unique(tab.Z_m_))==1, 'data not constant in z' );



% [x,ii] = sort(tab.X_m_)
% 
% u = tab.Velocity_i__m_s_(ii);
% 
% x = tab.X_m_(ii);
% 
% z = tab.Z_m_(ii);


siU = scatteredInterpolant(tab.X_m_,tab.Y_m_,tab.Velocity_i__m_s_);

figure, imagesc(tab.X_m_,tab.Y_m_,tab.Velocity_i__m_s_)