clear
clc

% deterministic parameters used by the struct "SP"
document = 'GSA_data_crop_V1.xlsx';

table_DP = readtable(document,'sheet','DP');

DP.capacity_pv = table_DP{ismember(table_DP{:,'variable'},'CAPACITY_PV'),'value'};
DP.capcost_ag = table_DP{ismember(table_DP{:,'variable'},'CAPCOST_AG'),'value'};
DP.prod_livestock = table_DP{ismember(table_DP{:,'variable'},'PROD_LIVESTOCK'),'value'};
DP.qty_co2 = table_DP{ismember(table_DP{:,'variable'},'QTY_CO2'),'value'};
DP.qty_water = table_DP{ismember(table_DP{:,'variable'},'QTY_WATER'),'value'};
DP.OM_labour = table_DP{ismember(table_DP{:,'variable'},'OM_LABOUR'),'value'};
DP.OM_ag = table_DP{ismember(table_DP{:,'variable'},'OM_AG'),'value'};
DP.OM_pv = table_DP{ismember(table_DP{:,'variable'},'OM_PV'),'value'};

% stochasitc parameters used by the struct "DP"
SRI_level = input('Input the Sunshine level (1/2/3)£º');
table_SP = readtable(document,'sheet','SP');

SP.capcost_pv      = table_SP{ismember(table_SP{:,'variable'},'CAPCOST_PV'),'baseline'};
SP.EFF_DECAY       = table_SP{ismember(table_SP{:,'variable'},'EFF_DECAY'),'baseline'};
SP.TIME_PV         = table_SP{ismember(table_SP{:,'variable'},['TIME_PV_',num2str(SRI_level)]),'baseline'};
SP.ele_pv          = DP.capacity_pv * ( 1 - SP.EFF_DECAY) * SP.TIME_PV * 1000; %MW to kW
SP.fit_pv          = table_SP{ismember(table_SP{:,'variable'},['FIT_PV_',num2str(SRI_level)]),'baseline'};
SP.profit_plant    = table_SP{ismember(table_SP{:,'variable'},'PROFIT_PLANT'),'baseline'};
SP.price_livestock = table_SP{ismember(table_SP{:,'variable'},'PRICE_LIVESTOCK'),'baseline'};
SP.qty_tour        = table_SP{ismember(table_SP{:,'variable'},'QTY_TOUR'),'baseline'}; 
SP.price_tour      = table_SP{ismember(table_SP{:,'variable'},'PRICE_TOUR'),'baseline'};
SP.price_co2       = table_SP{ismember(table_SP{:,'variable'},'PRICE_CO2'),'baseline'};
SP.price_water     = table_SP{ismember(table_SP{:,'variable'},'PRICE_WATER'),'baseline'};
SP.rent_land       = table_SP{ismember(table_SP{:,'variable'},'RENT_LAND'),'baseline'};
SP.tax_rate        = table_SP{ismember(table_SP{:,'variable'},'TAX_RATE'),'baseline'};
table_SP.Properties.RowNames = table_SP.variable;

salvage = 0.05 * SP.capcost_pv * 0.35;

netprofit = SP.ele_pv * SP.fit_pv ...
    - DP.qty_water * SP.price_water ...
    - SP.rent_land...
    - DP.OM_labour...
    - DP.OM_ag...
    - DP.OM_pv;

tax = netprofit * SP.tax_rate;

SPP = (SP.capcost_pv + DP.capcost_ag - salvage)/(SP.ele_pv * SP.fit_pv + ...
    SP.profit_plant + DP.prod_livestock * SP.price_livestock...
    + SP.qty_tour * SP.price_tour + DP.qty_co2 * SP.price_co2...
    - DP.qty_water * SP.price_water - SP.rent_land...
    - DP.OM_labour - DP.OM_ag - DP.OM_pv - tax);

if SRI_level == 1
    table_SP('TIME_PV_2',:) = [];
    table_SP('TIME_PV_3',:) = [];
    table_SP.Properties.RowNames{end} = 'TIME_PV';
    
    table_SP('FIT_PV_2',:) = [];
    table_SP('FIT_PV_3',:) = [];
    table_SP.Properties.RowNames{2} = 'FIT_PV';
    
elseif SRI_level == 2
    table_SP('TIME_PV_1',:) = [];
    table_SP('TIME_PV_3',:) = [];
    table_SP.Properties.RowNames{end} = 'TIME_PV';
    
    table_SP('FIT_PV_1',:) = [];
    table_SP('FIT_PV_3',:) = [];
    table_SP.Properties.RowNames{2} = 'FIT_PV';   
    
elseif SRI_level == 3
    table_SP('TIME_PV_1',:) = [];
    table_SP('TIME_PV_2',:) = [];
    table_SP.Properties.RowNames{end} = 'TIME_PV';
    
    table_SP('FIT_PV_1',:) = [];
    table_SP('FIT_PV_2',:) = [];
    table_SP.Properties.RowNames{2} = 'FIT_PV';       
end

disp('*****************************')
fprintf('The baseline of SPP is: %d \n', SPP);
disp('*****************************')

figure('name','SPP')

fig3 = figure('name','Sobol & Spearman');
fig4 = figure('name','Sort the indices');

alpha = 0.05;

for i = 1:4
    % uncertainty generation
    Perc = alpha * i;
    num_MC = 100000;

    Model.capcost_pv = makedist(table_SP{'CAPCOST_PV','distribution'}{1},...
        'mu',SP.capcost_pv ,...
        'sigma',Perc * SP.capcost_pv);
    Model.capcost_pv = truncate(Model.capcost_pv, ...
        table_SP{'CAPCOST_PV','upper'}, ...
        2 * table_SP{'CAPCOST_PV','upper'} - table_SP{'CAPCOST_PV','lower'});  
    Group.capcost_pv = random(Model.capcost_pv,[num_MC,1]);
    Group.capcost_pv = SP.capcost_pv - (Group.capcost_pv - SP.capcost_pv);

    Model.EFF_DECAY = makedist(table_SP{'EFF_DECAY','distribution'}{1},...
        'mu',SP.EFF_DECAY ,...
        'sigma',Perc * SP.EFF_DECAY);
    Model.EFF_DECAY = truncate(Model.EFF_DECAY, ...
        table_SP{'EFF_DECAY','lower'}, ...
        table_SP{'EFF_DECAY','upper'});
    Group.EFF_DECAY = random(Model.EFF_DECAY,[num_MC,1]);

    if SRI_level == 2 || SRI_level == 3
        Model.TIME_PV = makedist(table_SP{'TIME_PV','distribution'}{1},...
            'mu',SP.TIME_PV ,...
            'sigma',Perc * SP.TIME_PV);
        Model.TIME_PV = truncate(Model.TIME_PV, ...
            table_SP{'TIME_PV','lower'}, ...
            table_SP{'TIME_PV','upper'});      
        Group.TIME_PV = random(Model.TIME_PV, [num_MC,1]);
    elseif SRI_level == 1
        Model.TIME_PV = makedist(table_SP{'TIME_PV','distribution'}{1},...
            'mu',SP.TIME_PV ,...
            'sigma',Perc * SP.TIME_PV);
        Model.TIME_PV = truncate(Model.TIME_PV, ...
            table_SP{'TIME_PV','lower'}, ...
            table_SP{'TIME_PV','upper'});        
        Group.TIME_PV = random(Model.TIME_PV,[num_MC,1]);
    end   
    Group.ele_pv = DP.capacity_pv .* ( 1 - Group.EFF_DECAY) .* Group.TIME_PV * 1000; %MW to kW

    Model.fit_pv = makedist(table_SP{'FIT_PV','distribution'}{1},'mu',SP.fit_pv ,'sigma',Perc * SP.fit_pv);
    Model.fit_pv = truncate(Model.fit_pv, ...
        table_SP{'FIT_PV','upper'}, ...
        2 * table_SP{'FIT_PV','upper'} - table_SP{'FIT_PV','lower'});   
    Group.fit_pv = random(Model.fit_pv,[num_MC,1]);
    Group.fit_pv = SP.fit_pv - (Group.fit_pv - SP.fit_pv);

    Model.profit_plant = makedist(table_SP{'PROFIT_PLANT','distribution'}{1},...
        'mu',SP.profit_plant,...
        'sigma',Perc * SP.profit_plant);
    Model.profit_plant = truncate(Model.profit_plant, ...
        table_SP{'PROFIT_PLANT','upper'}, ...
        2 * table_SP{'PROFIT_PLANT','upper'} - table_SP{'PROFIT_PLANT','lower'});   
    Group.profit_plant = random(Model.profit_plant,[num_MC,1]);
    Group.profit_plant = SP.profit_plant - (Group.profit_plant - SP.profit_plant);
    
    Model.price_livestock = makedist(table_SP{'PRICE_LIVESTOCK','distribution'}{1},...
        'lower',table_SP{'PRICE_LIVESTOCK','lower'},...
        'upper',table_SP{'PRICE_LIVESTOCK','upper'});
    Group.price_livestock = random(Model.price_livestock,[num_MC,1]);

    Model.qty_tour = makedist(table_SP{'QTY_TOUR','distribution'}{1},...
        'mu',SP.qty_tour,...
        'sigma',Perc * SP.qty_tour);
    Model.qty_tour = truncate(Model.qty_tour, ...
        table_SP{'QTY_TOUR','lower'}, ...
        table_SP{'QTY_TOUR','upper'});   
    Group.qty_tour = random(Model.qty_tour,[num_MC,1]);

    Model.price_tour = makedist(table_SP{'PRICE_TOUR','distribution'}{1}, ...
        'mu',SP.price_tour, ...
        'sigma',Perc * SP.price_tour);
    Model.price_tour = truncate(Model.price_tour, ...
        table_SP{'PRICE_TOUR','lower'}, ...
        table_SP{'PRICE_TOUR','upper'});
    Group.price_tour = random(Model.price_tour,[num_MC,1]);
    
    Model.price_co2 = makedist(table_SP{'PRICE_CO2','distribution'}{1},...
        'mu',SP.price_co2,...
        'sigma',Perc * SP.price_co2);
    Model.price_co2 = truncate(Model.price_co2, ...
        table_SP{'PRICE_CO2','lower'}, ...
        table_SP{'PRICE_CO2','upper'});   
    Group.price_co2 = random(Model.price_co2, [num_MC,1]);
    
    Model.price_water = makedist(table_SP{'PRICE_WATER','distribution'}{1},...
        'mu',SP.price_water,...
        'sigma',Perc * SP.price_water);
    Model.price_water = truncate(Model.price_water, ...
        table_SP{'PRICE_WATER','lower'}, ...
        table_SP{'PRICE_WATER','upper'});      
    Group.price_water = random(Model.price_water,[num_MC,1]);

    Model.rent_land = makedist(table_SP{'RENT_LAND','distribution'}{1},...
        'mu',SP.rent_land,...
        'sigma',Perc * SP.rent_land);
    Model.rent_land = truncate(Model.rent_land, ...
        table_SP{'RENT_LAND','lower'}, ...
        table_SP{'RENT_LAND','upper'});    
    Group.rent_land = random(Model.rent_land,[num_MC,1]);

    Model.tax_rate = makedist(table_SP{'TAX_RATE','distribution'}{1},...
        'lower',table_SP{'TAX_RATE','lower'} ,...
        'upper',table_SP{'TAX_RATE','upper'});
    Model.tax_rate = truncate(Model.tax_rate, ...
        table_SP{'TAX_RATE','lower'}, ...
        table_SP{'TAX_RATE','upper'});    
    Group.tax_rate = random(Model.tax_rate,[num_MC,1]);

    Group.salvage = 0.05 * Group.capcost_pv * 0.35;

    Group.netprofit = Group.ele_pv .* Group.fit_pv ...
        - DP.qty_water .* Group.price_water ...
        - Group.rent_land...
        - DP.OM_labour...
        - DP.OM_ag...
        - DP.OM_pv;

    Group.tax = Group.netprofit .* Group.tax_rate;

    Group.SPP = (Group.capcost_pv + DP.capcost_ag - Group.salvage)./(Group.ele_pv .* Group.fit_pv + ...
        Group.profit_plant + DP.prod_livestock .* Group.price_livestock...
        + Group.qty_tour .* Group.price_tour + DP.qty_co2 .* Group.price_co2...
        - DP.qty_water .* Group.price_water - Group.rent_land...
        - DP.OM_labour - DP.OM_ag - DP.OM_pv - Group.tax);
    
    subplot(2,2,i)
    hist(Group.SPP)
%     title(['level: ', num2str(SRI_level)])
    title(num2str(Perc))
    hold on
    axis([0 10 -inf inf]) 
    drawnow
    
    disp('*****************************')
    fprintf('The mean of %d is: %d \n', [Perc, mean(Group.SPP)]);
    fprintf('The std of %d is: %d \n', [Perc, std(Group.SPP)]);
    disp('*****************************')

    % uncertainty parameter including: 
    % {capcost_pv, EFF_DECAY, TIME_PV, fit_pv, profit_plant, price_livestock, qty_tour
    % price_tour, price_co2, price_water, rent_land, tax_rate}
    
    GSA.matrix = [Group.capcost_pv, Group.EFF_DECAY, Group.TIME_PV,...
        Group.fit_pv, Group.profit_plant, Group.price_livestock,...
        Group.qty_tour, Group.price_tour, Group.price_co2,...
        Group.price_water, Group.rent_land, Group.tax_rate,...
        Group.SPP];

    GSA.coeff = corr(GSA.matrix, 'type', 'Spearman');

    figure('name','Coefficient')
    h = heatmap(GSA.coeff);

    Labels_list = {'capcost_pv', 'EFF_DECAY', 'TIME_PV',...
        'fit_pv', 'profit_plant', 'price_livestock',...
        'qty_tour', 'price_tour', 'price_co2',...
        'price_water', 'rent_land', 'tax_rate',...
        'SPP'};
    
    h.YDisplayLabels = Labels_list;
    h.XDisplayLabels = Labels_list;
    h.FontSize = 25;
    title(['level: ', num2str(SRI_level)])

    % Sobol_based Monte Carlo
    num_Sobol = num_MC / 2;
    
    Sobol_A.capcost_pv = Group.capcost_pv(1:num_Sobol);
    Sobol_B.capcost_pv = Group.capcost_pv(num_Sobol+1:num_MC);
    
    Sobol_A.price_livestock = Group.price_livestock(1:num_Sobol);
    Sobol_B.price_livestock = Group.price_livestock(num_Sobol+1:num_MC);
    
    Sobol_A.EFF_DECAY = Group.EFF_DECAY(1:num_Sobol);
    Sobol_B.EFF_DECAY = Group.EFF_DECAY(num_Sobol+1:num_MC);
 
    Sobol_A.TIME_PV = Group.TIME_PV(1:num_Sobol);
    Sobol_B.TIME_PV = Group.TIME_PV(num_Sobol+1:num_MC);
    
    Sobol_A.fit_pv = Group.fit_pv(1:num_Sobol);
    Sobol_B.fit_pv = Group.fit_pv(num_Sobol+1:num_MC);         
      
    Sobol_A.qty_tour = Group.qty_tour(1:num_Sobol);
    Sobol_B.qty_tour = Group.qty_tour(num_Sobol+1:num_MC);        
    
    Sobol_A.price_tour = Group.price_tour(1:num_Sobol);
    Sobol_B.price_tour = Group.price_tour(num_Sobol+1:num_MC);        

    Sobol_A.rent_land = Group.rent_land(1:num_Sobol);
    Sobol_B.rent_land = Group.rent_land(num_Sobol+1:num_MC);                
    
    Sobol_A.tax_rate = Group.tax_rate(1:num_Sobol);
    Sobol_B.tax_rate = Group.tax_rate(num_Sobol+1:num_MC);    
      
    Sobol_A.price_co2 = Group.price_co2(1:num_Sobol);
    Sobol_B.price_co2 = Group.price_co2(num_Sobol+1:num_MC);    
    
    Sobol_A.profit_plant = Group.profit_plant(1:num_Sobol);
    Sobol_B.profit_plant = Group.profit_plant(num_Sobol+1:num_MC);     
    
    Sobol_A.price_water = Group.price_water(1:num_Sobol);
    Sobol_B.price_water = Group.price_water(num_Sobol+1:num_MC);     
       
    Sobol_A = Crop_SPP_cal(Sobol_A, DP);
    Sobol_B = Crop_SPP_cal(Sobol_B, DP);
    
    A = [Sobol_A.capcost_pv, Sobol_A.EFF_DECAY, Sobol_A.TIME_PV, ...
        Sobol_A.fit_pv, Sobol_A.profit_plant, Sobol_A.price_livestock, ...
        Sobol_A.qty_tour, Sobol_A.price_tour, Sobol_A.price_co2, ...
        Sobol_A.price_water, Sobol_A.rent_land, Sobol_A.tax_rate];
    
    B = [Sobol_B.capcost_pv, Sobol_B.EFF_DECAY, Sobol_B.TIME_PV, ...
        Sobol_B.fit_pv, Sobol_B.profit_plant, Sobol_B.price_livestock, ...
        Sobol_B.qty_tour, Sobol_B.price_tour, Sobol_B.price_co2, ...
        Sobol_B.price_water, Sobol_B.rent_land, Sobol_B.tax_rate];
    
    [num_Sobol, dim] = size(A);
    
    AB = zeros(num_Sobol,dim,dim);
    for ii = 1:dim
        AB(:,:,ii) = A;
        AB(:,ii,ii) = B(:,ii);
    end
    
    Sobol_AB = struct([]);
    for ii = 1:dim
        Sobol_AB(ii).capcost_pv = AB(:,1,ii);
        Sobol_AB(ii).EFF_DECAY = AB(:,2,ii);
        Sobol_AB(ii).TIME_PV = AB(:,3,ii);
        Sobol_AB(ii).fit_pv = AB(:,4,ii);
        Sobol_AB(ii).profit_plant = AB(:,5,ii);  
        Sobol_AB(ii).price_livestock = AB(:,6,ii);
        Sobol_AB(ii).qty_tour = AB(:,7,ii);
        Sobol_AB(ii).price_tour = AB(:,8,ii);
        Sobol_AB(ii).price_co2 = AB(:,9,ii);
        Sobol_AB(ii).price_water = AB(:,10,ii);
        Sobol_AB(ii).rent_land = AB(:,11,ii);
        Sobol_AB(ii).tax_rate = AB(:,12,ii);  
    end
    
    for ii = 1:dim
        clear inter_struct
        inter_struct = Crop_SPP_cal(Sobol_AB(ii), DP);
        Sobol_AB(ii).SPP = inter_struct.SPP;
    end

    Var_Y = var(Sobol_A.SPP);
    % First order S1i
    S1 = zeros(1,dim);
    for ii = 1:dim
        sumvar = 0;
        for jj = 1:num_Sobol
            sumvar = sumvar + Sobol_B.SPP(jj)*(Sobol_AB(ii).SPP(jj)-Sobol_A.SPP(jj));
        end
        S1(1,ii) = sumvar/num_Sobol/Var_Y;
    end  
    % Total effect indices ST_i
            
    ST = zeros(1,dim);
    for ii = 1:dim
        sumvar = 0;
        for jj = 1:num_Sobol
            sumvar = sumvar + (Sobol_A.SPP(jj)-Sobol_AB(ii).SPP(jj))^2;
        end
        ST(1,ii) = sumvar/2/num_Sobol/Var_Y; 
    end
    
    sum_S1 = sum(S1);
    sum_ST = sum(ST);
    
    Labels = {'capcost_pv', 'EFF_DECAY', 'TIME_PV',...
        'fit_pv', 'profit_plant', 'price_livestock',...
        'qty_tour', 'price_tour', 'price_co2',...
        'price_water', 'rent_land', 'tax_rate'};
    
    figure(fig3) 
    subplot(2,2,1)
    plot(GSA.coeff(end,1:end-1))
    set(gca,'Xticklabel',Labels);
    title('Spearman')
    
    subplot(2,2,2)
    plot(abs(GSA.coeff(end,1:end-1)))
    set(gca,'Xticklabel',Labels);
    title('Absolute value of Spearman')
    
    subplot(2,2,3)
    plot(S1)
    set(gca,'Xticklabel',Labels);
    title('First order indices')
    
    subplot(2,2,4)
    plot(ST)
    set(gca,'Xticklabel',Labels);
    title('Total effect indices')
    
    Indices_table = table();
    Indices_table.variable = Labels';
    Indices_table.spearman = GSA.coeff(end,1:end-1)';
    Indices_table.abs_spearman = abs(GSA.coeff(end,1:end-1))';
    Indices_table.first_order = S1';
    Indices_table.total_effect = ST';
    
    disp('****************for spearman****************')
    sort_spearman = sortrows(Indices_table,'abs_spearman','descend');
    sort_spearman.rank = [1:height(sort_spearman)]';
    disp(sort_spearman(:,{'rank','variable','abs_spearman'}));
    
    disp('****************for first order sobol****************')
    sort_first_order = sortrows(Indices_table,'first_order','descend');
    sort_first_order.rank = [1:height(sort_first_order)]';
    disp(sort_first_order(:,{'rank','variable','first_order'})) 
    
    disp('****************for total effect sobol****************')
    sort_total_effect = sortrows(Indices_table,'total_effect','descend');
    sort_total_effect.rank = [1:height(sort_total_effect)]';
    disp(sort_total_effect(:,{'rank','variable','total_effect'}))
    
    figure(fig4) 
    subplot(3,1,1)
    plot(sort_spearman.abs_spearman) 
    set(gca,'Xticklabel',sort_spearman.variable,'FontSize',15);
    title('Absolute value of Spearman')
    
    subplot(3,1,2)
    plot(sort_first_order.first_order) 
    set(gca,'Xticklabel',sort_first_order.variable,'FontSize',15);
    title('First order indices')
    
    subplot(3,1,3)
    plot(sort_total_effect.total_effect)     
    set(gca,'Xticklabel',sort_total_effect.variable,'FontSize',15);
    title('Total effect indices')
      
    Perc = alpha * i;
     
    name_sheet = {'0p05','0p10','0p15','0p20'};
    name_excel_SPP = ['crop_SPP_level_', num2str(SRI_level), '.xlsx'];
    name_excel_sensit = ['crop_Sensit_level_', num2str(SRI_level), '.xlsx'];
    
    xlswrite(name_excel_SPP,Group.SPP,name_sheet{i});
    writetable(Indices_table,name_excel_sensit,'sheet',name_sheet{i});
    sum_sobol(i,:) = sum(Indices_table{:,{'first_order','total_effect'}});
end
close all
xlswrite(name_excel_SPP,SPP,'baseline');





