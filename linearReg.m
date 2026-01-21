function [PVx, PVy] = linearReg(curDetX, curDetY, curDetZ)
    curDetZ_=ones(size(curDetX));
    curDetZ_(:,1) = curDetZ_(:,1)*curDetZ(:,1);
    curDetZ_(:,2) = curDetZ_(:,2)*curDetZ(:,2);
    curDetZ_(:,3) = curDetZ_(:,3)*curDetZ(:,3);
    totalPions = size(curDetX,1);
    PVx = zeros(totalPions,1);
    PVy = zeros(totalPions,1);
    for i = 1:totalPions
        xyz = [curDetX(i,:);curDetY(i,:);curDetZ_(i,:)];
        
        xyz0 = mean(xyz,2);
        A = xyz-xyz0;
        [U,S,~] = svd(A);
        d = U(:,1);
        t = d'*A;
        t1 = min(t);
        t2 = max(t);
        xzyl = xyz0 + [t1,t2].*d; 
        
        t = -xzyl(3,1)./-xzyl(3,2);
        PVx(i) = xzyl(1,1) + t.*(xzyl(1,2)-xzyl(1,1));
        PVy(i) = xzyl(2,1) + t.*(xzyl(2,2)-xzyl(2,1));
%         % Check
%         x = xyz(1,:);
%         y = xyz(2,:);
%         z = xyz(3,:);
%         xl = xzyl(1,:);
%         yl = xzyl(2,:);
%         zl = xzyl(3,:);
% 
%         close all
%         hold on
%         plot3(x,y,z,'o');
%         plot3(xl,yl,zl,'r');
%         axis equal
    end
end