% Function Stokes_Continuity_solver()
% This function formulates and solves  
% Stokes and Continuity equations defined on 2D staggered irregularly spaced grid
% with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy)
% given distribution of right parts for all equations (RX,RY,RC) on the grid 
% and given variable shear (etas) and normal (etan) viscosity distributions 
% pressure is normalized relative to given value (prnorm) in the first cell
%
% Velocity Boundary condition specified by bleft,bright,btop,bbottom,bintern 
% are implemented from ghost nodes 
% directly into Stokes and continuity equations
%
% Function returns solution for velocity and pressure (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function[vx,resx,vy,resy,pr,resc]=Stokes_Continuity_solver(prfirst,etas,etan,xnum,ynum,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern)
% 
% Staggered Grid for Multigrid
% 
%     vx       vx       vx    
%
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Basic (density) nodes are shown with +
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Boundary conditions
% Pressure boundary condition
% Pressure in first cell
bpres=0;
prnorm=prfirst(2);
% Channel flow top->bottom
if (prfirst(1)==1)
    bpres=1;
    prnorm=prfirst(2);
end

% Velocity boundary conditions
btopx(:,1)=btop(:,1);
btopx(:,2)=btop(:,2);
btopy(:,1)=btop(:,3);
btopy(:,2)=btop(:,4);
bbottomx(:,1)=bbottom(:,1);
bbottomx(:,2)=bbottom(:,2);
bbottomy(:,1)=bbottom(:,3);
bbottomy(:,2)=bbottom(:,4);
bleftx(:,1)=bleft(:,1);
bleftx(:,2)=bleft(:,2);
blefty(:,1)=bleft(:,3);
blefty(:,2)=bleft(:,4);
brightx(:,1)=bright(:,1);
brightx(:,2)=bright(:,2);
brighty(:,1)=bright(:,3);
brighty(:,2)=bright(:,4);

% Prescribed internal velocity condition

% Computing grid steps for basic nodes
xstp=zeros(xnum-1,1);
ystp=zeros(ynum-1,1);
for i=1:1:xnum-1
    xstp(i)=gridx(i+1)-gridx(i);
end
for i=1:1:ynum-1
    ystp(i)=gridy(i+1)-gridy(i);
end
% Computing grid steps for vx and vy nodes
xstpc=zeros(xnum,1);
ystpc=zeros(ynum,1);
% First and last steps (for ghost nodes)
xstpc(1)=xstp(1);
xstpc(xnum)=xstp(xnum-1);
ystpc(1)=ystp(1);
ystpc(ynum)=ystp(ynum-1);
for i=2:1:xnum-1
    xstpc(i)=(gridx(i+1)-gridx(i-1))/2;
end
for i=2:1:ynum-1
    ystpc(i)=(gridy(i+1)-gridy(i-1))/2;
end

% Average x and y steps
xstpavr=(gridx(xnum)-gridx(1))/(xnum-1);
ystpavr=(gridy(ynum)-gridy(1))/(ynum-1);

% Koefficient for scaling pressure
pscale=2*etan(1)/(xstpavr+ystpavr);

% Horizontal shift index
ynum3=(ynum-1)*3;

% Creating matrix
L=sparse((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
R=zeros((xnum-1)*(ynum-1)*3,1);

% Solving of Stokes and continuity equations on nodes
for i=1:1:ynum-1
    for j=1:1:xnum-1
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        
        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum-1 && (j~=bintern(1) || i<bintern(2) || i>bintern(3)))
            % x-Stokes equation stensil
            %     +----------------------+----------------------+   
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i-1,j)                 ------------    
            %     |                      |                      |
            %     |                      |                      |
            %     +-----vy(i-1,j)---etas(i,j+1)---vy(i-1,j+1)----    ystpc(i)--------
            %     |                      |                      |
            %     |                      |                      |
            % vx(i,j-1)  pr(i,j)      vx(i,j)     P(i,j+1)   vx(i,j+1)-----   ystp(i)
            %     |     etan(i,j)        |       etan(i,j+1)    |
            %     |                      |                      |
            %     +------vy(i,j)---etas(i+1,j+1)---vy(i,j+1)-----    ystpc(i+1)------
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i+1,j)                 -----------    
            %     |                      |                      |
            %     |                      |                      |
            %     +----------------------+----------------------+   
            %     |         xstp(j)      |      xstp(j+1)       |   
            %                 |      xstpc(j+1)     |   
            % Right Part
            R(ivx,1)=RX(i+1,j+1);
            % Computing Current x-Stokes coefficients
            % Central Vx node
            L(ivx,ivx)=-2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1)-(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
            % Left Vx node
            if (j>1)
                ivxleft=ivx-ynum3;
                L(ivx,ivxleft)=2*etan(i,j)/xstp(j)/xstpc(j+1);
            else
                L(ivx,ivx)=L(ivx,ivx)+bleftx(i+1,2)*2*etan(i,j)/xstp(j)/xstpc(j+1);
                R(ivx,1)=R(ivx,1)-bleftx(i+1,1)*2*etan(i,j)/xstp(j)/xstpc(j+1);
            end
            % Right Vx node
            if (j<xnum-2)
                ivxright=ivx+ynum3;
                L(ivx,ivxright)=2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
            else
                L(ivx,ivx)=L(ivx,ivx)+brightx(i+1,2)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
                R(ivx,1)=R(ivx,1)-brightx(i+1,1)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
            end
            % Top Vx node
            if (i>1)
                ivxtop=ivx-3;
                L(ivx,ivxtop)=etas(i,j+1)/ystpc(i)/ystp(i);
            else
                L(ivx,ivx)=L(ivx,ivx)+btopx(j+1,2)*etas(i,j+1)/ystpc(i)/ystp(i);
                R(ivx,1)=R(ivx,1)-btopx(j+1,1)*etas(i,j+1)/ystpc(i)/ystp(i);
            end
            % Bottom Vx node
            if (i<ynum-1)
                ivxbottom=ivx+3;
                L(ivx,ivxbottom)=etas(i+1,j+1)/ystpc(i+1)/ystp(i);
            else
                L(ivx,ivx)=L(ivx,ivx)+bbottomx(j+1,2)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);
                R(ivx,1)=R(ivx,1)-bbottomx(j+1,1)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);
            end
            % Top Left Vy node
            if (i>1)
                ivytopleft=ivx-3+1;
                L(ivx,ivytopleft)=etas(i,j+1)/xstpc(j+1)/ystp(i);
            else
                ivybottomleft=ivx+1;
                L(ivx,ivybottomleft)=btopy(j+1,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
                R(ivx,1)=R(ivx,1)-btopy(j+1,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);
            end
            % Top Right Vy node
            if (i>1)
                ivytopright=ivx-3+1+ynum3;
                L(ivx,ivytopright)=-etas(i,j+1)/xstpc(j+1)/ystp(i);
            else
                ivybottomright=ivx+1+ynum3;
                L(ivx,ivybottomright)=-btopy(j+2,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
                R(ivx,1)=R(ivx,1)+btopy(j+2,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);
            end
            % Bottom Left Vy node
            if (i<ynum-1)
                ivybottomleft=ivx+1;
                if (i>1)
                    L(ivx,ivybottomleft)=-etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                else
                    L(ivx,ivybottomleft)=L(ivx,ivybottomleft)-etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                end
            else
                ivytopleft=ivx-3+1;
                L(ivx,ivytopleft)=L(ivx,ivytopleft)-bbottomy(j+1,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                R(ivx,1)=R(ivx,1)+bbottomy(j+1,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
            end
            % Bottom Right Vy node
            if (i<ynum-1)
                ivybottomright=ivx+1+ynum3;
                if (i>1)
                    L(ivx,ivybottomright)=etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                else
                    L(ivx,ivybottomright)=L(ivx,ivybottomright)+etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                end
            else
                ivytopright=ivx-3+1+ynum3;
                L(ivx,ivytopright)=L(ivx,ivytopright)+bbottomy(j+2,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
                R(ivx,1)=R(ivx,1)-bbottomy(j+2,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
            end
            % Left P node
            iprleft=ivx+2;
            L(ivx,iprleft)=pscale/xstpc(j+1);
            % Right P node
            iprright=ivx+2+ynum3;
            L(ivx,iprright)=-pscale/xstpc(j+1);
            
        % Ghost Vx_parameter=0 used for numbering, internal prescribed horizontal velocity
        else
            L(ivx,ivx)=2*pscale/(xstpavr+ystpavr);
            if (j~=bintern(1) || i<bintern(2) || i>bintern(3))
                R(ivx,1)=0;
            else
                % Internal prescribed horizontal velocity
                R(ivx,1)=2*pscale/(xstpavr+ystpavr)*bintern(4);
            end
                
        end 
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (i<ynum-1 && (j~=bintern(5) || i<bintern(6) || i>bintern(7)))
            % y-Stokes equation stensil
            %     +-------------------- -+-------vy(i-1,j)------+----------------------+-----    
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i,j-1)     P(i,j)    vx(i,j)                   |ystp(i)-------    
            %     |                      |        etan(i,j)     |                      |
            %     |                      |                      |                      |
            %     +-----vy(i,j-1)---etas(i+1,j)---vy(i,j)--etas(i+1,j+1)---vy(i,j+1)---+----- ystpc(i+1)
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i+1,j-1)  P(i+1,j)   vx(i+1,j)                 |ystp(i+1)------    
            %     |                      |      etan(i+1,j)     |                      |
            %     |                      |                      |                      |
            %     +----------------------+-------vy(i+1,j)------+----------------------+-----
            %               |          xstpc(j)      |         xstpc(j+1)      |   
            %                            |          xstp(j)     |
            %
            % Right Part
            R(ivy,1)=RY(i+1,j+1);
            % Computing Current y-Stokes coefficients
            % Central Vy node
            L(ivy,ivy)=-2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);
            % Top Vy node
            if(i>1)
                ivytop=ivy-3;
                L(ivy,ivytop)=2*etan(i,j)/ystp(i)/ystpc(i+1);
            else
                % Add boundary condition for the top Vy node
                L(ivy,ivy)=L(ivy,ivy)+btopy(j+1,2)*2*etan(i,j)/ystp(i)/ystpc(i+1);
                R(ivy,1)=R(ivy,1)-btopy(j+1,1)*2*etan(i,j)/ystp(i)/ystpc(i+1);
            end
            % Bottom Vy node
            if(i<ynum-2)
                ivybottom=ivy+3;
                L(ivy,ivybottom)=2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
            else
                % Add boundary condition for the bottom Vy node
                L(ivy,ivy)=L(ivy,ivy)+bbottomy(j+1,2)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
                R(ivy,1)=R(ivy,1)-bbottomy(j+1,1)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
            end
            % Left Vy node
            if(j>1)
                ivyleft=ivy-ynum3;
                L(ivy,ivyleft)=etas(i+1,j)/xstpc(j)/xstp(j);
            else
                % Add boundary condition for the left Vy node
                L(ivy,ivy)=L(ivy,ivy)+blefty(i+1,2)*etas(i+1,j)/xstpc(j)/xstp(j);
                R(ivy,1)=R(ivy,1)-blefty(i+1,1)*etas(i+1,j)/xstpc(j)/xstp(j);
            end
            % Right Vy node
            if(j<xnum-1)
                ivyright=ivy+ynum3;
                L(ivy,ivyright)=etas(i+1,j+1)/xstpc(j+1)/xstp(j);
            else
                % Add boundary condition for the right Vy node
                L(ivy,ivy)=L(ivy,ivy)+brighty(i+1,2)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);
                R(ivy,1)=R(ivy,1)-brighty(i+1,1)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);
            end
            % Top left Vx node
            if (j>1)
                ivxtopleft=ivy-1-ynum3;
                L(ivy,ivxtopleft)=etas(i+1,j)/ystpc(i+1)/xstp(j);
            else
                ivxtopright=ivy-1;
                L(ivy,ivxtopright)=bleftx(i+1,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
                R(ivy,1)=R(ivy,1)-bleftx(i+1,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);
            end
            % Bottom left Vx node
            if (j>1)
                ivxbottomleft=ivy-1+3-ynum3;
                L(ivy,ivxbottomleft)=-etas(i+1,j)/ystpc(i+1)/xstp(j);
            else
                ivxbottomright=ivy-1+3;
                L(ivy,ivxbottomright)=-bleftx(i+2,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
                R(ivy,1)=R(ivy,1)+bleftx(i+2,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);
            end
            % Top right Vx node
            if (j<xnum-1)
                ivxtopright=ivy-1;
                if(j>1)
                    L(ivy,ivxtopright)=-etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                else
                    L(ivy,ivxtopright)=L(ivy,ivxtopright)-etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                end
            else
                ivxtopleft=ivy-1-ynum3;
                L(ivy,ivxtopleft)=L(ivy,ivxtopleft)-brightx(i+1,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                R(ivy,1)=R(ivy,1)+brightx(i+1,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
           end
            % Bottom right Vx node
            if (j<xnum-1)
                ivxbottomright=ivy-1+3;
                if(j>1)
                    L(ivy,ivxbottomright)=etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                else
                    L(ivy,ivxbottomright)=L(ivy,ivxbottomright)+etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                end
            else
                ivxbottomleft=ivy-1+3-ynum3;
                L(ivy,ivxbottomleft)=L(ivy,ivxbottomleft)+brightx(i+2,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
                R(ivy,1)=R(ivy,1)-brightx(i+2,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
           end
            % Top P node
            iprtop=ivy+1;
            L(ivy,iprtop)=pscale/ystpc(i+1);
            % Bottom P node
            iprbottom=ivy+1+3;
            L(ivy,iprbottom)=-pscale/ystpc(i+1);
            
        % Ghost Vy_parameter=0 used for numbering
        else
            L(ivy,ivy)=2*pscale/(xstpavr+ystpavr);
            if (j~=bintern(5) || i<bintern(6) || i>bintern(7))
                R(ivy,1)=0;
            else
                % Internal prescribed horizontal velocity
                R(ivy,1)=2*pscale/(xstpavr+ystpavr)*bintern(8);
            end
        end

        
        % Continuity equation dvx/dx+dvy/dy=RC
        if ( ((j>1 || i>1) && bpres==0) || (i>1 && i<ynum-1 && bpres==1) || (j>1 && j<xnum-1 && bpres==2) ) 
            % Continuity equation stensil
            %     +-----vy(i-1,j)--------+--------
            %     |                      |
            %     |                      |
            % vx(i,j-1)  pr(i,j)      vx(i,j) ystp(i)
            %     |                      |
            %     |                      |
            %     +------vy(i,j)---------+--------
            %     |        xstp(j)       |
            %
            % Right Part
            R(ipr,1)=RC(i,j);
            % Computing Current Continuity coefficients
            % Left Vx node
            if (j>1)
                ivxleft=ipr-2-ynum3;
                L(ipr,ivxleft)=-pscale/xstp(j);
                % Add boundary condition for the right Vx node
                if (j==xnum-1)
                    L(ipr,ivxleft)=L(ipr,ivxleft)+brightx(i+1,2)*pscale/xstp(j);
                    R(ipr,1)=R(ipr,1)-brightx(i+1,1)*pscale/xstp(j);
                end
            end
            % Right Vx node
            if (j<xnum-1)
                ivxright=ipr-2;
                L(ipr,ivxright)=pscale/xstp(j);
                % Add boundary condition for the left Vx node
                if (j==1)
                    L(ipr,ivxright)=L(ipr,ivxright)-bleftx(i+1,2)*pscale/xstp(j);
                    R(ipr,1)=R(ipr,1)+bleftx(i+1,1)*pscale/xstp(j);
                end
            end
            % Top Vy node
            if (i>1)
                ivytop=ipr-1-3;
                L(ipr,ivytop)=-pscale/ystp(i);
                % Add boundary condition for the bottom Vy node
                if (i==ynum-1)
                    L(ipr,ivytop)=L(ipr,ivytop)+bbottomy(j+1,2)*pscale/ystp(i);
                    R(ipr,1)=R(ipr,1)-bbottomy(j+1,1)*pscale/ystp(i);
                end
            end
            % Bottom Vy node
            if (i<ynum-1)
                ivybottom=ipr-1;
                L(ipr,ivybottom)=pscale/ystp(i);
                % Add boundary condition for the top Vy node
                if (i==1)
                    L(ipr,ivybottom)=L(ipr,ivybottom)-btopy(j+1,2)*pscale/ystp(i);
                    R(ipr,1)=R(ipr,1)+btopy(j+1,1)*pscale/ystp(i);
                end
            end
            
        % Pressure definition for the boundary condition regions
        else
            % Pressure definition in one cell
            if (bpres==0)
                L(ipr,ipr)=2*pscale/(xstpavr+ystpavr);
                R(ipr,1)=2*prnorm/(xstpavr+ystpavr);
            end
            % Pressure definition at the top and bottom
            if (bpres==1)
                L(ipr,ipr)=2*pscale/(xstpavr+ystpavr);
                if (i==1)
                    R(ipr,1)=2*prnorm/(xstpavr+ystpavr);
                else
                    R(ipr,1)=0;
                end
            end
            % Pressure definition at the left and right
            if (bpres==2)
                L(ipr,ipr)=2*pscale/(xstpavr+ystpavr);
                if (j==1)
                    R(ipr,1)=2*prnorm/(xstpavr+ystpavr);
                else
                    R(ipr,1)=0;
                end
            end
        end
             
    end            
end

S=L\R;

% Reload solution
vx=zeros(ynum+1,xnum);
vy=zeros(ynum,xnum+1);
pr=zeros(ynum-1,xnum-1);
for i=1:1:ynum-1
    for j=1:1:xnum-1
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        % Reload Vx
        if (j<xnum-1)
            vx(i+1,j+1)=S(ivx);
        end
        % Reload Vy
        if (i<ynum-1)
            vy(i+1,j+1)=S(ivy);
        end
        % Reload P
        pr(i,j)=S(ipr)*pscale;
    end
end

% Apply vx boundary conditions
% Left,Right Boundary
for i=1:1:ynum+1
    vx(i,1)=bleftx(i,1)+bleftx(i,2)*vx(i,2);
    vx(i,xnum)=brightx(i,1)+brightx(i,2)*vx(i,xnum-1);
end
% Top, Bottom Boundary
for j=1:1:xnum
    vx(1,j)=btopx(j,1)+btopx(j,2)*vx(2,j);
    vx(ynum+1,j)=bbottomx(j,1)+bbottomx(j,2)*vx(ynum,j);
end

% Apply vy boundary conditions
% Left,Right Boundary
for i=1:1:ynum
    vy(i,1)=blefty(i,1)+blefty(i,2)*vy(i,2);
    vy(i,xnum+1)=brighty(i,1)+brighty(i,2)*vy(i,xnum);
end
% Top, Bottom Boundary
for j=1:1:xnum+1
    vy(1,j)=btopy(j,1)+btopy(j,2)*vy(2,j);
    vy(ynum,j)=bbottomy(j,1)+bbottomy(j,2)*vy(ynum-1,j);
end

% Initialize residual arrays
resx=zeros(ynum+1,xnum);
resy=zeros(ynum,xnum+1);
resc=zeros(ynum-1,xnum-1);
% Computing residuals
for i=1:1:ynum+1;
    for j=1:1:xnum+1;

        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum+1 && (j~=bintern(1) || i<bintern(2) || i>bintern(3)))
            % vx-Boundrary conditions 
            if (i==1 || i==ynum+1 || j==1 || j==xnum);
                resx(i,j)=0;
            % x-Stokes equation
            % x-Stokes equation stensil
            %     +----------------------+----------------------+   
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i-1,j)                 ------------    
            %     |                      |                      |
            %     |                      |                      |
            %     +-----vy(i-1,j)---etas(i-1,j)---vy(i-1,j+1)----    ystpc(i-1)----
            %     |                      |                      |
            %     |                      |                      |
            % vx(i,j-1)  pr(i-1,j-1)  vx(i,j)     P(i-1,j)   vx(i,j+1)-----   ystp(i-1)
            %     |    etan(i-1,j-1)     |       etan(i-1,j)    |
            %     |                      |                      |
            %     +-------vy(i,j)----etas(i,j)-----vy(i,j+1)-----    ystpc(i)------
            %     |                      |                      |
            %     |                      |                      |
            %     |                   vx(i+1,j)                 -----------    
            %     |                      |                      |
            %     |                      |                      |
            %     +----------------------+----------------------+   
            %     |         xstp(j-1)    |      xstp(j)         |   
            %                 |      xstpc(j)        |   
            else
                % Computing Current x-Stokes residual
                % dSIGMAxx/dx-dP/dx
                resx(i,j)=RX(i,j)-(2*(etan(i-1,j)*(vx(i,j+1)-vx(i,j))/xstp(j)-etan(i-1,j-1)*(vx(i,j)-vx(i,j-1))/xstp(j-1))-(pr(i-1,j)-pr(i-1,j-1)))/xstpc(j);
                % dSIGMAxy/dy
                resx(i,j)=resx(i,j)-(etas(i,j)*((vx(i+1,j)-vx(i,j))/ystpc(i)+(vy(i,j+1)-vy(i,j))/xstpc(j))-etas(i-1,j)*((vx(i,j)-vx(i-1,j))/ystpc(i-1)+(vy(i-1,j+1)-vy(i-1,j))/xstpc(j)))/ystp(i-1);
            end
        end
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (i<ynum+1 && (j~=bintern(5) || i<bintern(6) || i>bintern(7)))
            % vy-Boundrary conditions 
            if (i==1 || i==ynum || j==1 || j==xnum+1);
                resy(i,j)=0;
            %y-Stokes equation
            % y-Stokes equation stensil
            %     +-------------------- -+-------vy(i-1,j)------+----------------------+-----    
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i,j-1)  P(i-1,j-1)   vx(i,j)                   |ystp(i-1)-------    
            %     |                      |      etan(i-1,j-1)   |                      |
            %     |                      |                      |                      |
            %     +-----vy(i,j-1)---etas(i,j-1)----vy(i,j)----etas(i,j)---vy(i,j+1)----+----- ystpc(i)
            %     |                      |                      |                      |
            %     |                      |                      |                      |
            %     |                  vx(i+1,j-1)  P(i,j-1)   vx(i+1,j)                 |ystp(i)------    
            %     |                      |      etan(i,j-1)     |                      |
            %     |                      |                      |                      |
            %     +----------------------+-------vy(i+1,j)------+----------------------+-----
            %               |        xstpc(j-1)      |        xstpc(j)      |   
            %                            |       xstp(j-1)      |
            %
            else
                % Computing current residual
                % dSIGMAyy/dy-dP/dy
                resy(i,j)=RY(i,j)-(2*(etan(i,j-1)*(vy(i+1,j)-vy(i,j))/ystp(i)-etan(i-1,j-1)*(vy(i,j)-vy(i-1,j))/ystp(i-1))-(pr(i,j-1)-pr(i-1,j-1)))/ystpc(i);
                % dSIGMAxy/dx
                resy(i,j)=resy(i,j)-(etas(i,j)*((vy(i,j+1)-vy(i,j))/xstpc(j)+(vx(i+1,j)-vx(i,j))/ystpc(i))-etas(i,j-1)*((vy(i,j)-vy(i,j-1))/xstpc(j-1)+(vx(i+1,j-1)-vx(i,j-1))/ystpc(i)))/xstp(j-1);
            end
        end
            
        % Continuity equation dvx/dx+dvy/dy=RC
        if (i<ynum && j<xnum);
            % Continuity equation stensil
            %     +------vy(i,j+1)-------+--------
            %     |                      |
            %     |                      |
            % vx(i+1,j)   pr(i,j)   vx(i+1,j+1) ystp(i)
            %     |                      |
            %     |                      |
            %     +-----vy(i+1,j+1)------+--------
            %     |        xstp(j)       |
            %
            % Computing current residual
            resc(i,j)=RC(i,j)-((vx(i+1,j+1)-vx(i+1,j))/xstp(j)+(vy(i+1,j+1)-vy(i,j+1))/ystp(i));
        end
             
    end            
end


