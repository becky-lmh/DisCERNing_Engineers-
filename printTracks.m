function printTracks(xPos,yPos,zPos,printing)

if printing == true
    plot3(xPos',yPos',zPos');
    xlabel("x (m)")
    ylabel("y (m)")
    zlabel("z (m)")
    grid on
    set(gca,'FontSize',10)
end

end