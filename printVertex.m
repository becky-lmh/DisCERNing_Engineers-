function printVertex(pionCharge, printing)

[~,~,ic] = unique(table2array(pionCharge(:,"collisionNum")));
pionPerVertex = accumarray(ic,1);

[C,~,ic] = unique(pionPerVertex(:,1));
a_counts = accumarray(ic,1);
if printing == true 
    figure()
    bar(C,a_counts);
    grid on
    xlabel("Charged pions produced per vertex")
    ylabel("Frequency")
    set(gca,'FontSize',10)
    %title("Frequency distribution of number of charged pions per vertex")
end
end