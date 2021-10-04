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
function[vx,resx,vy,resy,pr,resc]=Stokes_Continuity_solver_ll(prfirst,etas,etan,xnum,ynum,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern)
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


% % Creating matrix
% L=sparse((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
% R=zeros((xnum-1)*(ynum-1)*3,1);

Lcx1 = cell(ynum,xnum); Lcx2 = cell(ynum,xnum);
Lcx3 = cell(ynum,xnum); Lcx4 = cell(ynum,xnum);
Lcx5 = cell(ynum,xnum); Lcx6 = cell(ynum,xnum);
Lcx7 = cell(ynum,xnum); Lcx8 = cell(ynum,xnum);
Lcx9 = cell(ynum,xnum); Lcx10 = cell(ynum,xnum);
Lcx11 = cell(ynum,xnum);
Rcx = cell(ynum,xnum);

Lcy1 = cell(ynum,xnum); Lcy2 = cell(ynum,xnum);
Lcy3 = cell(ynum,xnum); Lcy4 = cell(ynum,xnum);
Lcy5 = cell(ynum,xnum); Lcy6 = cell(ynum,xnum);
Lcy7 = cell(ynum,xnum); Lcy8 = cell(ynum,xnum);
Lcy9 = cell(ynum,xnum); Lcy10 = cell(ynum,xnum);
Lcy11 = cell(ynum,xnum);
Rcy = cell(ynum,xnum);

Lc1 = cell(ynum,xnum); Lc2 = cell(ynum,xnum);
Lc3 = cell(ynum,xnum); Lc4 = cell(ynum,xnum);
Lc5 = cell(ynum,xnum); Lc6 = cell(ynum,xnum);
Lc7 = cell(ynum,xnum); Lc8 = cell(ynum,xnum);
Lc9 = cell(ynum,xnum); Lc10 = cell(ynum,xnum);
Lc11 = cell(ynum,xnum);
Rc = cell(ynum,xnum);

temp_i_limit = ynum-1;
temp_j_limit = xnum-1;


% Solving of Stokes and continuity equations on nodes
parfor i=1:temp_i_limit
    temp_i = i; % Se lasci i come indice (l'indice del parfor) lui se lo 
                % perde. Devi rispecificarlo dopo il parfor affinché ogni 
                % worker (o core) abbia ben chiaro l'indice i su cui sta
                % lavorando, in questo caso temp_i.
    for j=1:temp_j_limit
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(temp_i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        
        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum-1 && (j~=bintern(1) || temp_i<bintern(2) || temp_i>bintern(3)))
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
            
            % Gli indici della cella di output DEVONO essere i e j. Poi
            % all'interno della cella crei un vettore con 3 elementi: la
            % coordinata riga (qui ivx), la coordinata colonna (qui 1) e il
            % valore della variabile (qui RX(temp_i+1,j+1))
            
            % Right Part
            Rcx{i,j} = [ivx,1,RX(temp_i+1,j+1)];
            % Computing Current x-Stokes coefficients
            % Central Vx node
            Lcx1{i,j} = [ivx,ivx,-2*(etan(temp_i,j+1)/xstp(j+1)+etan(temp_i,j)/xstp(j))/xstpc(j+1)-(etas(temp_i+1,j+1)/ystpc(temp_i+1)+etas(temp_i,j+1)/ystpc(temp_i))/ystp(temp_i)];
            % Left Vx node
            if (j>1)
                ivxleft=ivx-ynum3;
                Lcx2{i,j} = [ivx,ivxleft,2*etan(temp_i,j)/xstp(j)/xstpc(j+1)];
            else
                temp_Lcx1 = Lcx1{i,j}(3);
                Lcx1{i,j} = [ivx,ivx,temp_Lcx1+bleftx(temp_i+1,2)*2*etan(temp_i,j)/xstp(j)/xstpc(j+1)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-bleftx(temp_i+1,1)*2*etan(temp_i,j)/xstp(j)/xstpc(j+1)];
            end
            % Right Vx node
            if (j<xnum-2)
                ivxright=ivx+ynum3;
                Lcx3{i,j} = [ivx,ivxright,2*etan(temp_i,j+1)/xstp(j+1)/xstpc(j+1)];
            else
                temp_Lcx1 = Lcx1{i,j}(3);
                Lcx1{i,j} = [ivx,ivx,temp_Lcx1+brightx(temp_i+1,2)*2*etan(temp_i,j+1)/xstp(j+1)/xstpc(j+1)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-brightx(temp_i+1,1)*2*etan(temp_i,j+1)/xstp(j+1)/xstpc(j+1)];
            end
            % Top Vx node
            if (temp_i>1)
                ivxtop=ivx-3;
                Lcx4{i,j} = [ivx,ivxtop,etas(temp_i,j+1)/ystpc(temp_i)/ystp(temp_i)];
            else
                temp_Lcx1 = Lcx1{i,j}(3);
                Lcx1{i,j} = [ivx,ivx,temp_Lcx1+btopx(j+1,2)*etas(temp_i,j+1)/ystpc(temp_i)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-btopx(j+1,1)*etas(temp_i,j+1)/ystpc(temp_i)/ystp(temp_i)];
            end
            % Bottom Vx node
            if (temp_i<ynum-1)
                ivxbottom=ivx+3;
                Lcx5{i,j} = [ivx,ivxbottom,etas(temp_i+1,j+1)/ystpc(temp_i+1)/ystp(temp_i)];
            else
                temp_Lcx1 = Lcx1{i,j}(3);
                Lcx1{i,j} = [ivx,ivx,temp_Lcx1+bbottomx(j+1,2)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-bbottomx(j+1,1)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/ystp(temp_i)];
            end
            % Top Left Vy node
            if (temp_i>1)
                ivytopleft=ivx-3+1;
                Lcx6{i,j} = [ivx,ivytopleft,etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
            else
                ivybottomleft=ivx+1;
                Lcx7{i,j} = [ivx,ivybottomleft,btopy(j+1,2)*etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-btopy(j+1,1)*etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
            end
            % Top Right Vy node
            if (temp_i>1)
                ivytopright=ivx-3+1+ynum3;
                Lcx8{i,j} = [ivx,ivytopright,-etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
            else
                ivybottomright=ivx+1+ynum3;
                Lcx9{i,j} = [ivx,ivybottomright,-btopy(j+2,2)*etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx+btopy(j+2,1)*etas(temp_i,j+1)/xstpc(j+1)/ystp(temp_i)];
            end
            % Bottom Left Vy node
            if (temp_i<ynum-1)
                ivybottomleft=ivx+1;
                if (temp_i>1)
                    Lcx7{i,j} = [ivx,ivybottomleft,-etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                else
                    temp_Lcx7 = Lcx7{i,j}(3);
                    Lcx7{i,j} = [ivx,ivybottomleft,temp_Lcx7-etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                end
            else
                ivytopleft=ivx-3+1;
                temp_Lcx6 = Lcx6{i,j}(3);
                Lcx6{i,j} = [ivx,ivytopleft,temp_Lcx6-bbottomy(j+1,2)*etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx+bbottomy(j+1,1)*etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
            end
            % Bottom Right Vy node
            if (temp_i<ynum-1)
                ivybottomright=ivx+1+ynum3;
                if (temp_i>1)
                    Lcx9{i,j} = [ivx,ivybottomright,etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                else
                    temp_Lcx9 = Lcx9{i,j}(3);
                    Lcx9{i,j} = [ivx,ivybottomright,temp_Lcx9+etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                end
            else
                ivytopright=ivx-3+1+ynum3;
                temp_Lcx8 = Lcx8{i,j}(3);
                Lcx8{i,j} = [ivx,ivytopright,temp_Lcx8+bbottomy(j+2,2)*etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
                temp_Rcx = Rcx{i,j}(3);
                Rcx{i,j} = [ivx,1,temp_Rcx-bbottomy(j+2,1)*etas(temp_i+1,j+1)/xstpc(j+1)/ystp(temp_i)];
            end
            % Left P node
            iprleft=ivx+2;
            Lcx10{i,j} = [ivx,iprleft,pscale/xstpc(j+1)];
            % Right P node
            iprright=ivx+2+ynum3;
            Lcx11{i,j} = [ivx,iprright,-pscale/xstpc(j+1)];
            
        % Ghost Vx_parameter=0 used for numbering, internal prescribed horizontal velocity
        else
            Lcx1{i,j} = [ivx,ivx,2*pscale/(xstpavr+ystpavr)];
            if (j~=bintern(1) || temp_i<bintern(2) || temp_i>bintern(3))
                Rcx{i,j} = [ivx,1,0];
            else
                % Internal prescribed horizontal velocity
                Rcx{i,j} = [ivx,1,2*pscale/(xstpavr+ystpavr)*bintern(4)];
            end
                
        end

            
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (temp_i<ynum-1 && (j~=bintern(5) || temp_i<bintern(6) || temp_i>bintern(7)))
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
            Rcy{i,j} = [ivy,1,RY(temp_i+1,j+1)];
            % Computing Current y-Stokes coefficients
            % Central Vy node
            Lcy1{i,j} = [ivy,ivy,-2*(etan(temp_i+1,j)/ystp(temp_i+1)+etan(temp_i,j)/ystp(temp_i))/ystpc(temp_i+1)-(etas(temp_i+1,j+1)/xstpc(j+1)+etas(temp_i+1,j)/xstpc(j))/xstp(j)];
            % Top Vy node
            if(temp_i>1)
                ivytop=ivy-3;
                Lcy2{i,j} = [ivy,ivytop,2*etan(temp_i,j)/ystp(temp_i)/ystpc(temp_i+1)];
            else
                % Add boundary condition for the top Vy node
                temp_Lcy1 = Lcy1{i,j}(3);
                Lcy1{i,j} = [ivy,ivy,temp_Lcy1+btopy(j+1,2)*2*etan(temp_i,j)/ystp(temp_i)/ystpc(temp_i+1)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-btopy(j+1,1)*2*etan(temp_i,j)/ystp(temp_i)/ystpc(temp_i+1)];
            end
            % Bottom Vy node
            if(temp_i<ynum-2)
                ivybottom=ivy+3;
                Lcy3{i,j} = [ivy,ivybottom,2*etan(temp_i+1,j)/ystp(temp_i+1)/ystpc(temp_i+1)];
            else
                % Add boundary condition for the bottom Vy node
                temp_Lcy1 = Lcy1{i,j}(3);
                Lcy1{i,j} = [ivy,ivy,temp_Lcy1+bbottomy(j+1,2)*2*etan(temp_i+1,j)/ystp(temp_i+1)/ystpc(temp_i+1)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-bbottomy(j+1,1)*2*etan(temp_i+1,j)/ystp(temp_i+1)/ystpc(temp_i+1)];
            end
            % Left Vy node
            if(j>1)
                ivyleft=ivy-ynum3;
                Lcy4{i,j} = [ivy,ivyleft,etas(temp_i+1,j)/xstpc(j)/xstp(j)];
            else
                % Add boundary condition for the left Vy node
                temp_Lcy1 = Lcy1{i,j}(3);
                Lcy1{i,j} = [ivy,ivy,temp_Lcy1+blefty(temp_i+1,2)*etas(temp_i+1,j)/xstpc(j)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-blefty(temp_i+1,1)*etas(temp_i+1,j)/xstpc(j)/xstp(j)];
            end
            % Right Vy node
            if(j<xnum-1)
                ivyright=ivy+ynum3;
                Lcy5{i,j} = [ivy,ivyright,etas(temp_i+1,j+1)/xstpc(j+1)/xstp(j)];
            else
                % Add boundary condition for the right Vy node
                temp_Lcy1 = Lcy1{i,j}(3);
                Lcy1{i,j} = [ivy,ivy,temp_Lcy1+brighty(temp_i+1,2)*etas(temp_i+1,j+1)/xstpc(j+1)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-brighty(temp_i+1,1)*etas(temp_i+1,j+1)/xstpc(j+1)/xstp(j)];
            end
            % Top left Vx node
            if (j>1)
                ivxtopleft=ivy-1-ynum3;
                Lcy6{i,j} = [ivy,ivxtopleft,etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
            else
                ivxtopright=ivy-1;
                Lcy7{i,j} = [ivy,ivxtopright,bleftx(temp_i+1,2)*etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-bleftx(temp_i+1,1)*etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
            end
            % Bottom left Vx node
            if (j>1)
                ivxbottomleft=ivy-1+3-ynum3;
                Lcy8{i,j} = [ivy,ivxbottomleft,-etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
            else
                ivxbottomright=ivy-1+3;
                Lcy9{i,j} = [ivy,ivxbottomright,-bleftx(temp_i+2,2)*etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy+bleftx(temp_i+2,1)*etas(temp_i+1,j)/ystpc(temp_i+1)/xstp(j)];
            end
            % Top right Vx node
            if (j<xnum-1)
                ivxtopright=ivy-1;
                if(j>1)
                    Lcy7{i,j} = [ivy,ivxtopright,-etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                else
                    temp_Lcy7 = Lcy7{i,j}(3);
                    Lcy7{i,j} = [ivy,ivxtopright,temp_Lcy7-etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                end
            else
                ivxtopleft=ivy-1-ynum3;
                temp_Lcy6 = Lcy6{i,j}(3);
                Lcy6{i,j} = [ivy,ivxtopleft,temp_Lcy6-brightx(temp_i+1,2)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy+brightx(temp_i+1,1)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
           end
            % Bottom right Vx node
            if (j<xnum-1)
                ivxbottomright=ivy-1+3;
                if(j>1)
                    Lcy9{i,j} = [ivy,ivxbottomright,etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                else
                    temp_Lcy9 = Lcy9{i,j}(3);
                    Lcy9{i,j} = [ivy,ivxbottomright,temp_Lcy9+etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                end
            else
                ivxbottomleft=ivy-1+3-ynum3;
                temp_Lcy8 = Lcy8{i,j}(3);
                Lcy8{i,j} = [ivy,ivxbottomleft,temp_Lcy8+brightx(temp_i+2,2)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
                temp_Rcy = Rcy{i,j}(3);
                Rcy{i,j} = [ivy,1,temp_Rcy-brightx(temp_i+2,1)*etas(temp_i+1,j+1)/ystpc(temp_i+1)/xstp(j)];
           end
            % Top P node
            iprtop=ivy+1;
            Lcy10{i,j} = [ivy,iprtop,pscale/ystpc(temp_i+1)];
            % Bottom P node
            iprbottom=ivy+1+3;
            Lcy11{i,j} = [ivy,iprbottom,-pscale/ystpc(temp_i+1)];
            
        % Ghost Vy_parameter=0 used for numbering
        else
            Lcy1{i,j} = [ivy,ivy,2*pscale/(xstpavr+ystpavr)];
            if (j~=bintern(5) || temp_i<bintern(6) || temp_i>bintern(7))
                Rcy{i,j} = [ivy,1,0];
            else
                % Internal prescribed horizontal velocity
                Rcy{i,j} = [ivy,1,2*pscale/(xstpavr+ystpavr)*bintern(8)];
            end
        end

        
        % Continuity equation dvx/dx+dvy/dy=RC
        if ( ((j>1 || temp_i>1) && bpres==0) || (temp_i>1 && temp_i<ynum-1 && bpres==1) || (j>1 && j<xnum-1 && bpres==2) ) 
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
            Rc{i,j} = [ipr,1,RC(temp_i,j)];
            % Computing Current Continuity coefficients
            % Left Vx node
            if (j>1)
                ivxleft=ipr-2-ynum3;
                Lc1{i,j} = [ipr,ivxleft,-pscale/xstp(j)];
                % Add boundary condition for the right Vx node
                if (j==xnum-1)
                    temp_Lc1 = Lc1{i,j}(3);
                    Lc1{i,j} = [ipr,ivxleft,temp_Lc1+brightx(temp_i+1,2)*pscale/xstp(j)];
                    temp_Rc = Rc{i,j}(3);
                    Rc{i,j} = [ipr,1,temp_Rc-brightx(temp_i+1,1)*pscale/xstp(j)];
                end
            end
            % Right Vx node
            if (j<xnum-1)
                ivxright=ipr-2;
                Lc2{i,j} = [ipr,ivxright,pscale/xstp(j)];
                % Add boundary condition for the left Vx node
                if (j==1)
                    temp_Lc2 = Lc2{i,j}(3);
                    Lc2{i,j} = [ipr,ivxright,temp_Lc2-bleftx(temp_i+1,2)*pscale/xstp(j)];
                    temp_Rc = Rc{i,j}(3);
                    Rc{i,j} = [ipr,1,temp_Rc+bleftx(temp_i+1,1)*pscale/xstp(j)];
                end
            end
            % Top Vy node
            if (temp_i>1)
                ivytop=ipr-1-3;
                Lc3{i,j} = [ipr,ivytop,-pscale/ystp(temp_i)];
                % Add boundary condition for the bottom Vy node
                if (temp_i==ynum-1)
                    temp_Lc3 = Lc3{i,j}(3);
                    Lc3{i,j} = [ipr,ivytop,temp_Lc3+bbottomy(j+1,2)*pscale/ystp(temp_i)];
                    temp_Rc = Rc{i,j}(3);
                    Rc{i,j} = [ipr,1,temp_Rc-bbottomy(j+1,1)*pscale/ystp(temp_i)];
                end
            end
            % Bottom Vy node
            if (temp_i<ynum-1)
                ivybottom=ipr-1;
                Lc4{i,j} = [ipr,ivybottom,pscale/ystp(temp_i)];
                % Add boundary condition for the top Vy node
                if (temp_i==1)
                    temp_Lc4 = Lc4{i,j}(3);
                    Lc4{i,j} = [ipr,ivybottom,temp_Lc4-btopy(j+1,2)*pscale/ystp(temp_i)];
                    temp_Rc = Rc{i,j}(3);
                    Rc{i,j} = [ipr,1,temp_Rc+btopy(j+1,1)*pscale/ystp(temp_i)];
                end
            end
            
        % Pressure definition for the boundary condition regions
        else
            % Pressure definition in one cell
            if (bpres==0)
                Lc5{i,j} = [ipr,ipr,2*pscale/(xstpavr+ystpavr)];
                Rc{i,j} = [ipr,1,2*prnorm/(xstpavr+ystpavr)];
            end
            % Pressure definition at the top and bottom
            if (bpres==1)
                Lc5{i,j} = [ipr,ipr,2*pscale/(xstpavr+ystpavr)];
                if (temp_i==1)
                    Rc{i,j} = [ipr,1,2*prnorm/(xstpavr+ystpavr)];
                else
                    Rc{i,j} = [ipr,1,0];
                end
            end
            % Pressure definition at the left and right
            if (bpres==2)
                Lc5{i,j} = [ipr,ipr,2*pscale/(xstpavr+ystpavr)];
                if (j==1)
                    Rc{i,j} = [ipr,1,2*prnorm/(xstpavr+ystpavr)];
                else
                    Rc{i,j} = [ipr,1,0];
                end
            end
        end
             
    end            
end

Lcx1 = Lcx1(~cellfun(@isempty,Lcx1)); Lcx2 = Lcx2(~cellfun(@isempty,Lcx2));
Lcx3 = Lcx3(~cellfun(@isempty,Lcx3)); Lcx4 = Lcx4(~cellfun(@isempty,Lcx4));
Lcx5 = Lcx5(~cellfun(@isempty,Lcx5)); Lcx6 = Lcx6(~cellfun(@isempty,Lcx6));
Lcx7 = Lcx7(~cellfun(@isempty,Lcx7)); Lcx8 = Lcx8(~cellfun(@isempty,Lcx8));
Lcx9 = Lcx9(~cellfun(@isempty,Lcx9)); Lcx10 = Lcx10(~cellfun(@isempty,Lcx10));
Lcx11 = Lcx11(~cellfun(@isempty,Lcx11));
Lcx1 = cell2mat(Lcx1); Lcx2 = cell2mat(Lcx2); Lcx3 = cell2mat(Lcx3);
Lcx4 = cell2mat(Lcx4); Lcx5 = cell2mat(Lcx5); Lcx6 = cell2mat(Lcx6);
Lcx7 = cell2mat(Lcx7); Lcx8 = cell2mat(Lcx8); Lcx9 = cell2mat(Lcx9);
Lcx10 = cell2mat(Lcx10); Lcx11 = cell2mat(Lcx11);
Lcx = [Lcx1;Lcx2;Lcx3;Lcx4;Lcx5;Lcx6;Lcx7;Lcx8;Lcx9;Lcx10;Lcx11];

Lcy1 = Lcy1(~cellfun(@isempty,Lcy1)); Lcy2 = Lcy2(~cellfun(@isempty,Lcy2));
Lcy3 = Lcy3(~cellfun(@isempty,Lcy3)); Lcy4 = Lcy4(~cellfun(@isempty,Lcy4));
Lcy5 = Lcy5(~cellfun(@isempty,Lcy5)); Lcy6 = Lcy6(~cellfun(@isempty,Lcy6));
Lcy7 = Lcy7(~cellfun(@isempty,Lcy7)); Lcy8 = Lcy8(~cellfun(@isempty,Lcy8));
Lcy9 = Lcy9(~cellfun(@isempty,Lcy9)); Lcy10 = Lcy10(~cellfun(@isempty,Lcy10));
Lcy11 = Lcy11(~cellfun(@isempty,Lcy11));
Lcy1 = cell2mat(Lcy1); Lcy2 = cell2mat(Lcy2); Lcy3 = cell2mat(Lcy3);
Lcy4 = cell2mat(Lcy4); Lcy5 = cell2mat(Lcy5); Lcy6 = cell2mat(Lcy6);
Lcy7 = cell2mat(Lcy7); Lcy8 = cell2mat(Lcy8); Lcy9 = cell2mat(Lcy9);
Lcy10 = cell2mat(Lcy10); Lcy11 = cell2mat(Lcy11);
Lcy = [Lcy1;Lcy2;Lcy3;Lcy4;Lcy5;Lcy6;Lcy7;Lcy8;Lcy9;Lcy10;Lcy11];

Lc1 = Lc1(~cellfun(@isempty,Lc1)); Lc2 = Lc2(~cellfun(@isempty,Lc2));
Lc3 = Lc3(~cellfun(@isempty,Lc3)); Lc4 = Lc4(~cellfun(@isempty,Lc4));
Lc5 = Lc5(~cellfun(@isempty,Lc5)); Lc6 = Lc6(~cellfun(@isempty,Lc6));
Lc7 = Lc7(~cellfun(@isempty,Lc7)); Lc8 = Lc8(~cellfun(@isempty,Lc8));
Lc9 = Lc9(~cellfun(@isempty,Lc9)); Lc10 = Lc10(~cellfun(@isempty,Lc10));
Lc11 = Lc11(~cellfun(@isempty,Lc11));
Lc1 = cell2mat(Lc1); Lc2 = cell2mat(Lc2); Lc3 = cell2mat(Lc3);
Lc4 = cell2mat(Lc4); Lc5 = cell2mat(Lc5); Lc6 = cell2mat(Lc6);
Lc7 = cell2mat(Lc7); Lc8 = cell2mat(Lc8); Lc9 = cell2mat(Lc9);
Lc10 = cell2mat(Lc10); Lc11 = cell2mat(Lc11);
Lc = [Lc1;Lc2;Lc3;Lc4;Lc5;Lc6;Lc7;Lc8;Lc9;Lc10;Lc11];

LL = [Lcx;Lcy;Lc];
L = sparse(LL(:,1),LL(:,2),LL(:,3),(xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
Rc = [Rcx(:);Rcy(:);Rc(:)]; Rc = cell2mat(Rc); Rc = sortrows(Rc);
R = Rc(:,end);

% Solve matrix
S=L\R;

% SparS = S; LparS = L; RparS = R;
% save('S_par_Stokes','SparS'); save('L_par_Stokes','LparS'); save('R_par_Stokes','RparS');

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
resx=zeros(ynum+1,xnum+1); temp_resx = zeros(ynum+1,xnum+1);
resy=zeros(ynum+1,xnum+1); temp_resy = zeros(ynum+1,xnum+1);
resc=zeros(ynum+1,xnum+1);
% Computing residuals
temp_i_limit = ynum+1;
temp_j_limit = xnum+1;
parfor i=1:temp_i_limit
    temp_i = i;
    for j=1:temp_j_limit
        % x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
        if (j<xnum+1 && (j~=bintern(1) || temp_i<bintern(2) || temp_i>bintern(3)))
            % vx-Boundrary conditions 
            if (temp_i==1 || temp_i==ynum+1 || j==1 || j==xnum)
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
                temp_resx(i,j)=RX(temp_i,j)-(2*(etan(temp_i-1,j)*(vx(temp_i,j+1)-vx(temp_i,j))/xstp(j)-etan(temp_i-1,j-1)*(vx(temp_i,j)-vx(temp_i,j-1))/xstp(j-1))-(pr(temp_i-1,j)-pr(temp_i-1,j-1)))/xstpc(j);
                % dSIGMAxy/dy
                resx(i,j)=temp_resx(i,j)-(etas(temp_i,j)*((vx(temp_i+1,j)-vx(temp_i,j))/ystpc(temp_i)+(vy(temp_i,j+1)-vy(temp_i,j))/xstpc(j))-etas(temp_i-1,j)*((vx(temp_i,j)-vx(temp_i-1,j))/ystpc(temp_i-1)+(vy(temp_i-1,j+1)-vy(temp_i-1,j))/xstpc(j)))/ystp(temp_i-1);
            end
        end
            
        % y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
        if (temp_i<ynum+1 && (j~=bintern(5) || temp_i<bintern(6) || temp_i>bintern(7)))
            % vy-Boundrary conditions 
            if (temp_i==1 || temp_i==ynum || j==1 || j==xnum+1)
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
                temp_resy(i,j)=RY(temp_i,j)-(2*(etan(temp_i,j-1)*(vy(temp_i+1,j)-vy(temp_i,j))/ystp(temp_i)-etan(temp_i-1,j-1)*(vy(temp_i,j)-vy(temp_i-1,j))/ystp(temp_i-1))-(pr(temp_i,j-1)-pr(temp_i-1,j-1)))/ystpc(temp_i);
                % dSIGMAxy/dx
                resy(i,j)=temp_resy(i,j)-(etas(temp_i,j)*((vy(temp_i,j+1)-vy(temp_i,j))/xstpc(j)+(vx(temp_i+1,j)-vx(temp_i,j))/ystpc(temp_i))-etas(temp_i,j-1)*((vy(temp_i,j)-vy(temp_i,j-1))/xstpc(j-1)+(vx(temp_i+1,j-1)-vx(temp_i,j-1))/ystpc(temp_i)))/xstp(j-1);
            end
        end
            
        % Continuity equation dvx/dx+dvy/dy=RC
        if (temp_i<ynum && j<xnum)
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
            resc(i,j)=RC(temp_i,j)-((vx(temp_i+1,j+1)-vx(temp_i+1,j))/xstp(j)+(vy(temp_i+1,j+1)-vy(temp_i,j+1))/ystp(temp_i));
        end
             
    end            
end

% ResxparS = resx; ResyparS = resy; RescparS = resy;
% save('Resx_par_Stokes','ResxparS'); save('Resy_par_Stokes','ResyparS'); save('Resc_par_Stokes','RescparS')