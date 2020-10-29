clear all
close all

%% Run multiple simulations to get average activations
iter=100;
num=0;
for i=1:iter    
    num=num+1;
    disp(num);
    
    results();
    
    filename=('PPNactivationslast.mat'); 
    m = matfile(filename,'Writable',true);
    new_row_of_values = [lc_tot_activ, lc_stn_tot_activ, lc_stn_on_tot_activ, lc_snr_tot_activ, lc_snr_on_tot_activ];
    if isprop(m, 'PPNactivationslast')
        s = size(m, 'PPNactivationslast');
        m.PPNactivationslast(s(1)+1, :) = new_row_of_values;
    else
        m.PPNactivationslast = new_row_of_values;
    end 
end

