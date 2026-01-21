function printNumPions(newZPos,numPions, inc, printing)

if printing == true
    figure()
    plot(newZPos(1,:),numPions)
    grid on
    xlabel("Distance from center of collision /m")
    ylabel("Number of pions")
    set(gca,'FontSize',10)
    %title("Pions detected along the beam pipe axis")
    
    figure()
    grad = (numPions(:,2:end) - numPions(:,1:end-1))./inc;
    grad = smoothdata(grad,'gaussian');
    plot(newZPos(1,2:end-1),grad(:,2:end))
    grid on
    xlabel("Distance from center of collison /m")
    ylabel("Gradient")
    %title("Gradient of pions detected graph")
% 
%     figure()
%     grad2 = (grad(:,2:end) - grad(:,1:end-1))./inc;
%     grad2 = smoothdata(grad2,'gaussian');
%     plot(newZPos(1,2:end-2),grad2(:,2:end))
%     grid on
%     xlabel("Distance from center of collison /m")
%     ylabel("Second Gradient")
    %title("Second Gradient of pions detected graph")
end

end