% Function Temperature_solver_grid()
% This function formulates and solves  
% Heat conservation equation defined on 2D irregularly spaced grid
% with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy)
% given distribution of right parts for all equations (RT) on the grid 
% and given RHO*CP and k values
%
% Thermal Boundary condition specified by arrays bleft(),bright(),btop(),bbot() 
%
% Function returns solution for new temperature (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[tknew,rest]=Temperature_solver_grid(timestep,xnum,ynum,gridx,gridy,kt,rhocp,tk,RT,bleft,bright,btop,bbottom)
% 
% Staggered Grid for Multigrid
%
%     T--------T--------T 
%     |        |        |
%     |        |        |
%     |        |        |
%     T--------T--------T 
%     |        |        |
%     |        |        |
%     |        |        |
%     T--------T--------T 
% 
% Lines show basic grid


% Computing grid steps for basic nodes
xstp=zeros(xnum-1,1);
ystp=zeros(ynum-1,1);
for i=1:1:xnum-1
    xstp(i)=gridx(i+1)-gridx(i);
end
for i=1:1:ynum-1
    ystp(i)=gridy(i+1)-gridy(i);
end
% Computing grid steps for qx and qy nodes
xstpc=zeros(xnum-2,1);
ystpc=zeros(ynum-2,1);
for i=1:1:xnum-2
    xstpc(i)=(gridx(i+2)-gridx(i))/2;
end
for i=1:1:ynum-2
    ystpc(i)=(gridy(i+2)-gridy(i))/2;
end


% Creating matrix
L=sparse(xnum*ynum,xnum*ynum);
R=zeros(xnum*ynum,1);

% Solving of Stokes and continuity equations on nodes
for i=1:1:ynum
    for j=1:1:xnum
        % Index for T
        itk=(j-1)*ynum+i;
        
        % Boundary conditions
        if (j==1 || j==xnum || i==1 || i==ynum)
            % Upper boundary: tk(i,j)=btop(1)+btop(2)*tk(i+1,j)
            if (i==1 && j>1 && j<xnum)
                % Right part
                R(itk,1)=btop(j,1);
                % Left part: 1*tk(i,j)-btop(2)*tk(i+1,j)
                L(itk,itk)=1;
                L(itk,itk+1)=-btop(j,2);
            end
            % Lower boundary: tk(i,j)=bbottom(1)+bbottom(2)*tk(i-1,j)
            if (i==ynum && j>1 && j<xnum)
                % Right part
                R(itk,1)=bbottom(j,1);
                % Left part: 1*tk(i,j)-bbottom(2)*tk(i-1,j)
                L(itk,itk)=1;
                L(itk,itk-1)=-bbottom(j,2);
            end
            % Left boundary: tk(i,j)=bleft(1)+bleft(2)*tk(i,j+1)
            if (j==1)
                % Right part
                R(itk,1)=bleft(i,1);
                % Left part: 1*tk(i,j)-bleft(2)*tk(i,j+1)
                L(itk,itk)=1;
                L(itk,itk+ynum)=-bleft(i,2);
            end
            % Right boundary: tk(i,j)=bright(1)+bright(2)*tk(i,j-1)
            if (j==xnum)
                % Right part
                R(itk,1)=bright(i,1);
                % Left part: 1*tk(i,j)-bright(2)*tk(i,j-1)
                L(itk,itk)=1;
                L(itk,itk-ynum)=-bright(i,2);
            end
            
        % Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
        else
            % Temperature equation stensil
            %
            %           | xstpc(j-1) |  
            %     | xstp(j-1) |  xstp(j)  |
            %     +-------tk(i-1,j)-------+--------- 
            %     |       kt(i-1,j)       |
            %     |           |           |
            %     |           |           |    ystp(i-1)-----
            %     |           |           |
            %     |           |           |
            % tk(i,j-1)----tk(i,j)-----tk(i,j+1)------  ystpc(i-1) 
            % kt(i,j-1)    kt(i,j)     kt(i,j+1) 
            %     |       rhocp(i,j)      |
            %     |           |           |    ystp(i)-------
            %     |           |           |
            %     |           |           |
            %     +-------tk(i+1,j)-------+---------
            %             kt(i+1,j)
            %
            % Right Part
            R(itk,1)=RT(i,j)+tk(i,j)*rhocp(i,j)/timestep;
            % Computing coefficients for the left part
            % Central T node
            L(itk,itk)=rhocp(i,j)/timestep+((kt(i,j-1)+kt(i,j))/xstp(j-1)+(kt(i,j)+kt(i,j+1))/xstp(j))/2/xstpc(j-1)+((kt(i-1,j)+kt(i,j))/ystp(i-1)+(kt(i,j)+kt(i+1,j))/ystp(i))/2/ystpc(i-1);
            % Left T node
            L(itk,itk-ynum)=-(kt(i,j-1)+kt(i,j))/2/xstp(j-1)/xstpc(j-1);
            % Right T node
            L(itk,itk+ynum)=-(kt(i,j)+kt(i,j+1))/2/xstp(j)/xstpc(j-1);
            % Upper T node
            L(itk,itk-1)=-(kt(i-1,j)+kt(i,j))/2/ystp(i-1)/ystpc(i-1);
            % Lower T node
            L(itk,itk+1)=-(kt(i,j)+kt(i+1,j))/2/ystp(i)/ystpc(i-1);
        end
             
    end            
end


% Solve matrix
S=L\R;

% Reload solution
for i=1:1:ynum
    for j=1:1:xnum
        % Index for T
        itk=(j-1)*ynum+i;
        % Reload T
        tknew(i,j)=S(itk);
    end
end


% Computing residuals
for i=1:1:ynum;
    for j=1:1:xnum;

        % Boundary conditions
        if (j==1 || j==xnum || i==1 || i==ynum)
            rest(i,j)=0;
        % Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
        else
            % Computing Current Temperature equation residual
            % Ht-DT/dt
            rest(i,j)=RT(i,j)-rhocp(i,j)*(tknew(i,j)-tk(i,j))/timestep;
            % -dqx/dx
            rest(i,j)=rest(i,j)+((kt(i,j)+kt(i,j+1))*(tknew(i,j+1)-tknew(i,j))/xstp(j)-(kt(i,j-1)+kt(i,j))*(tknew(i,j)-tknew(i,j-1))/xstp(j-1))/2/xstpc(j-1);
            % -dqy/dy
            rest(i,j)=rest(i,j)+((kt(i,j)+kt(i+1,j))*(tknew(i+1,j)-tknew(i,j))/ystp(i)-(kt(i-1,j)+kt(i,j))*(tknew(i,j)-tknew(i-1,j))/ystp(i-1))/2/ystpc(i-1);
        end
    end
end


