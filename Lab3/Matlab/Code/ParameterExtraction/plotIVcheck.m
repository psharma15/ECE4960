% ------------------------------------------------------------------------
% Visual validation. Id vs Vds curve
% Takes input as the paramters and variables and sMeasured
% ------------------------------------------------------------------------
% Pragya Sharma, March 28 2017
% ------------------------------------------------------------------------
function plotIVcheck(sol,sMeas)
    %% plot sMeasured
    startRow = 1;
    endRow = 101;
    totRow = 1010;
    sMod = ekvModel(sol,sMeas);
    figure
    while endRow <= totRow
        tempColor = rand(1,3);
        plot(sMeas(startRow:endRow,2),sMod(startRow:endRow),'color',tempColor);
        hold on
        plot(sMeas(startRow:endRow,2),sMeas(startRow:endRow,3),'--','color',tempColor);
        text(sMeas(endRow,2),sMod(endRow),['V_g_s = ', num2str(sMeas(endRow,1))],'FontSize',10);
        hold on
        startRow = endRow + 1;
        endRow = endRow + 101;
    end
    
    legend('S_m_o_d_e_l','S_m_e_a_s_u_r_e_d');
    axis ([0 6 0 6e-3])
    ylabel('I_d','FontSize',12)
    xlabel('V_d_s','FontSize',12)
    title('Plot of I_d vs V_d_s','FontSize',12)
    
    %% Plot sModel with different Vds and Vgs
    Vds = [1;4];
    m = length(Vds);
    Vgs = (0.1:0.1:3)';
    n = length(Vgs);
    Vgs = repmat(Vgs,m,1);
    Vds = repmat(Vds,1,n)';
    Vds = Vds(:);
    measData = [Vgs,Vds];
    sMod = ekvModel(sol,measData);
    startRow = 1;
    endRow = n;
    totRow = n*m;
    figure
    while endRow <= totRow
        tempColor = rand(1,3);
        semilogy(measData(startRow:endRow,1),sMod(startRow:endRow),'color',tempColor);
        text(measData(endRow,1),sMod(endRow),['V_d_s = ', num2str(measData(endRow,2))],'FontSize',10);
        hold on
        startRow = endRow + 1;
        endRow = endRow + n;
    end
    
    axis ([0 3.5 0 1e-2])
    ylabel('log(I_d)','FontSize',12)
    xlabel('V_g_s','FontSize',12)
    title('Plot of log(I_d) vs V_g_s','FontSize',12)
    
    %% Check and Validation
    % --------------------------------------------------------------------
    tol = 1e-2;
    Vt = 26e-3;
    % --------------------------------------------------------------------
    % Vgs < Vth, Vds > 3*Vth: Turn-Off region
    startPt = 31;
    endPt = 35;
    sMod1 = sMod(startPt:endPt);
    idTurnOff = zeros(endPt-startPt+1,1);
    for i = 1:endPt-startPt+1
        idTurnOff(i) = sol(1)*exp((sol(2)*(measData(i,1)-sol(3)))/Vt)*(1 - exp(-measData(i,2)/Vt));
    end
    chkId1 = norm(idTurnOff - sMod1);
    if chkId1 < tol
        fprintf('Solution converging to exponential function in Turn-off region.\n');
    end
    
    % --------------------------------------------------------------------
    % Vgs > Vth, Vds < Vds_saturation: Above threshold Linear
    startPt = 25;
    endPt = 30;
    sMod1 = sMod(startPt:endPt);
    idAboveThLin = zeros(endPt-startPt+1,1);
    for i = 1:endPt-startPt+1
        idAboveThLin(i) = (sol(1)/(4*Vt^2))*((sol(2)*(measData(i,1)-sol(3)))^2 - (sol(2)*(measData(i,1)-sol(3)) - measData(i,2))^2);
    end
    chkId2 = norm(idAboveThLin - sMod1);
    if chkId2 < tol
        fprintf('Solution converging to quadratic function in Above Threshold Linear region.\n');
    end

    % --------------------------------------------------------------------
    % Vgs > Vth, Vds > Vdsat: Above Threshold Saturation
    startPt = 44;
    endPt = 50;
    sMod1 = sMod(startPt:endPt);
    idAboveThSat = zeros(endPt-startPt+1,1);
    for i = 1:endPt-startPt+1
        idAboveThSat(i) = (sol(1)/(4*Vt^2))*(sol(2)*(measData(i,1)-sol(3)))^2;
    end
    chkId2 = norm(idAboveThSat - sMod1);
    if chkId2 < tol
        fprintf('Solution converging to quadratic function in Above Threshold Saturation region.\n');
    end

end