function plotSecondary(collision_tab,printing)

[C,~,ic] = unique(table2array(collision_tab(:,"Secondary")));
a_counts = accumarray(ic,1);
if printing == true
    figure()
    bar(C,a_counts);
    grid on
    xlabel("Secondary particle")
    ylabel("Number of particles produced")
    set(gca,'FontSize',10)
    %title("Distribution of secondary particles produced in 7 TeV gas-beam collision")
    %set(gca, 'yscale','log');
end

end