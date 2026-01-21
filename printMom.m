function printMom(pionCharge,printing)

if printing == true
    histogram(table2array(pionCharge(:,6)),20);
    set(gca,'yscale','log')
    grid on
    xlabel("Momentum / MeV")
    ylabel("Number of charged pions")
    %title("Momentum histogram for charged pions")
end

end