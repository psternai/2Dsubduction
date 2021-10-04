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
function[tknew,rest]=Temperature_solver_grid_ll(timestep,xnum,ynum,gridx,gridy,kt,rhocp,tk,RT,bleft,bright,btop,bbottom)
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


% % Creating matrix
% L=sparse(xnum*ynum,xnum*ynum);
% R=zeros(xnum*ynum,1);

Lc1 = cell(ynum,xnum); Lc2 = cell(ynum,xnum);
Lc3 = cell(ynum,xnum); Lc4 = cell(ynum,xnum);
Lc5 = cell(ynum,xnum); Lc6 = cell(ynum,xnum);
Lc7 = cell(ynum,xnum);
Rc = cell(ynum,xnum);

% Solving of Stokes and continuity equations on nodes
for i=1:ynum  %ynum è il numero di righe della matrice
    temp_i = i; % Se lasci i come indice (l'indice del parfor) lui se lo 
                % perde. Devi rispecificarlo dopo il parfor affinché ogni 
                % worker (o core) abbia ben chiaro l'indice i su cui sta
                % lavorando, in questo caso temp_i.
    for j=1:xnum
        % Index for T
        itk=(j-1)*ynum+temp_i;
%         cnt = cnt+1;
        % Boundary conditions
        if (j==1 || j==xnum || temp_i==1 || temp_i==ynum)
            % Upper boundary: tk(i,j)=btop(1)+btop(2)*tk(i+1,j)
            if (temp_i==1 && j>1 && j<xnum)
                % Right part
                Rc{i,j} = [itk,1,btop(j,1)];
                % Left part: 1*tk(i,j)-btop(2)*tk(i+1,j)
                Lc1{i,j}=[itk,itk,1];
%                 cnt = cnt+1
                Lc2{i,j}=[itk,itk+1,-btop(j,2)];
            end
            % Lower boundary: tk(i,j)=bbottom(1)+bbottom(2)*tk(i-1,j)
            if (temp_i==ynum && j>1 && j<xnum)
                % Right part
                Rc{i,j} = [itk,1,bbottom(j,1)];
                % Left part: 1*tk(i,j)-bbottom(2)*tk(i-1,j)
                Lc1{i,j} = [itk,itk,1];
%                 cnt = cnt+1
                Lc2{i,j} = [itk,itk-1,-bbottom(j,2)];
            end
            % Left boundary: tk(i,j)=bleft(1)+bleft(2)*tk(i,j+1)
            if (j==1)
                % Right part
                Rc{i,j} = [itk,1,bleft(temp_i,1)];
                % Left part: 1*tk(i,j)-bleft(2)*tk(i,j+1)
                Lc1{i,j} = [itk,itk,1];
%                 cnt = cnt+1
                Lc2{i,j} = [itk,itk+ynum,-bleft(temp_i,2)];
            end
            % Right boundary: tk(i,j)=bright(1)+bright(2)*tk(i,j-1)
            if (j==xnum)
                % Right part
                Rc{i,j} = [itk,1,bright(temp_i,1)];
                % Left part: 1*tk(i,j)-bright(2)*tk(i,j-1)
                Lc1{i,j} = [itk,itk,1];
%                 cnt = cnt+1
                Lc2{i,j} = [itk,itk-ynum,-bright(temp_i,2)];
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
            Rc{i,j} = [itk,1,RT(temp_i,j)+tk(temp_i,j)*rhocp(temp_i,j)/timestep];
            % Computing coefficients for the left part
            % Central T node
            Lc3{i,j} = [itk,itk,rhocp(temp_i,j)/timestep+((kt(temp_i,j-1)+kt(temp_i,j))/xstp(j-1)+(kt(temp_i,j)+kt(temp_i,j+1))/xstp(j))/2/xstpc(j-1)+((kt(temp_i-1,j)+kt(temp_i,j))/ystp(temp_i-1)+(kt(temp_i,j)+kt(temp_i+1,j))/ystp(temp_i))/2/ystpc(temp_i-1)];
%             cnt = cnt+1
            % Left T node
            Lc4{i,j} = [itk,itk-ynum,-(kt(temp_i,j-1)+kt(temp_i,j))/2/xstp(j-1)/xstpc(j-1)];
%             cnt = cnt+1
            % Right T node
            Lc5{i,j} = [itk,itk+ynum,-(kt(temp_i,j)+kt(temp_i,j+1))/2/xstp(j)/xstpc(j-1)];
%             cnt = cnt+1
            % Upper T node
            Lc6{i,j} = [itk,itk-1,-(kt(temp_i-1,j)+kt(temp_i,j))/2/ystp(temp_i-1)/ystpc(temp_i-1)];
%             cnt = cnt+1
            % Lower T node
            Lc7{i,j} = [itk,itk+1,-(kt(temp_i,j)+kt(temp_i+1,j))/2/ystp(temp_i)/ystpc(temp_i-1)];
        end
             
    end
%     Ll=sparse(xnum*ynum,xnum*ynum);
%     Rr=zeros(xnum*ynum,1);
%     temp_i = i; % Se lasci i come indice (l'indice del parfor) lui se lo 
%                 % perde. Devi rispecificarlo dopo il parfor affinché ogni 
%                 % worker (o core) abbia ben chiaro l'indice i su cui sta
%                 % lavorando, in questo caso temp_i.
%     for j=1:xnum
%         % Index for T
%         itk=(j-1)*ynum+temp_i;
%         % Boundary conditions
%         if (j==1 || j==xnum || temp_i==1 || temp_i==ynum)
%             % Upper boundary: tk(i,j)=btop(1)+btop(2)*tk(i+1,j)
%             if (temp_i==1 && j>1 && j<xnum)
%                 % Right part
%                 Rr(itk,1)=btop(j,1);
%                 % Left part: 1*tk(i,j)-btop(2)*tk(i+1,j)
%                 Ll(itk,itk)=1;
%                 Ll(itk,itk+1)=-btop(j,2);
%             end
%             % Lower boundary: tk(i,j)=bbottom(1)+bbottom(2)*tk(i-1,j)
%             if (temp_i==ynum && j>1 && j<xnum)
%                 % Right part
%                 Rr(itk,1)=bbottom(j,1);
%                 % Left part: 1*tk(i,j)-bbottom(2)*tk(i-1,j)
%                 Ll(itk,itk)=1;
%                 Ll(itk,itk-1)=-bbottom(j,2);
%             end
%             % Left boundary: tk(i,j)=bleft(1)+bleft(2)*tk(i,j+1)
%             if (j==1)
%                 % Right part
%                 Rr(itk,1)=bleft(temp_i,1);
%                 % Left part: 1*tk(i,j)-bleft(2)*tk(i,j+1)
%                 Ll(itk,itk)=1;
%                 Ll(itk,itk+ynum)=-bleft(temp_i,2);
%             end
%             % Right boundary: tk(i,j)=bright(1)+bright(2)*tk(i,j-1)
%             if (j==xnum)
%                 % Right part
%                 Rr(itk,1)=bright(temp_i,1);
%                 % Left part: 1*tk(i,j)-bright(2)*tk(i,j-1)
%                 Ll(itk,itk)=1;
%                 Ll(itk,itk-ynum)=-bright(temp_i,2);
%             end
%             
%         % Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
%         else
%             % Temperature equation stensil
%             %
%             %           | xstpc(j-1) |  
%             %     | xstp(j-1) |  xstp(j)  |
%             %     +-------tk(i-1,j)-------+--------- 
%             %     |       kt(i-1,j)       |
%             %     |           |           |
%             %     |           |           |    ystp(i-1)-----
%             %     |           |           |
%             %     |           |           |
%             % tk(i,j-1)----tk(i,j)-----tk(i,j+1)------  ystpc(i-1) 
%             % kt(i,j-1)    kt(i,j)     kt(i,j+1) 
%             %     |       rhocp(i,j)      |
%             %     |           |           |    ystp(i)-------
%             %     |           |           |
%             %     |           |           |
%             %     +-------tk(i+1,j)-------+---------
%             %             kt(i+1,j)
%             %
%             % Right Part
%             Rr(itk,1)=RT(temp_i,j)+tk(temp_i,j)*rhocp(temp_i,j)/timestep;
%             % Computing coefficients for the left part
%             % Central T node
%             Ll(itk,itk)=rhocp(temp_i,j)/timestep+((kt(temp_i,j-1)+kt(temp_i,j))/xstp(j-1)+(kt(temp_i,j)+kt(temp_i,j+1))/xstp(j))/2/xstpc(j-1)+((kt(temp_i-1,j)+kt(temp_i,j))/ystp(temp_i-1)+(kt(temp_i,j)+kt(temp_i+1,j))/ystp(temp_i))/2/ystpc(temp_i-1);
%             % Left T node
%             Ll(itk,itk-ynum)=-(kt(temp_i,j-1)+kt(temp_i,j))/2/xstp(j-1)/xstpc(j-1);
%             % Right T node
%             Ll(itk,itk+ynum)=-(kt(temp_i,j)+kt(temp_i,j+1))/2/xstp(j)/xstpc(j-1);
%             % Upper T node
%             Ll(itk,itk-1)=-(kt(temp_i-1,j)+kt(temp_i,j))/2/ystp(temp_i-1)/ystpc(temp_i-1);
%             % Lower T node
%             Ll(itk,itk+1)=-(kt(temp_i,j)+kt(temp_i+1,j))/2/ystp(temp_i)/ystpc(temp_i-1);
%             Lc{i,j} = Ll;
%             Rc{i,j} = Rr;
%         end
%              
%     end
end

Lc1 = Lc1(~cellfun(@isempty,Lc1)); Lc2 = Lc2(~cellfun(@isempty,Lc2));
Lc3 = Lc3(~cellfun(@isempty,Lc3)); Lc4 = Lc4(~cellfun(@isempty,Lc4));
Lc5 = Lc5(~cellfun(@isempty,Lc5)); Lc6 = Lc6(~cellfun(@isempty,Lc6));
Lc7 = Lc7(~cellfun(@isempty,Lc7));
Lc1 = cell2mat(Lc1); Lc2 = cell2mat(Lc2); Lc3 = cell2mat(Lc3);
Lc4 = cell2mat(Lc4); Lc5 = cell2mat(Lc5); Lc6 = cell2mat(Lc6);
Lc7 = cell2mat(Lc7);
Lc = [Lc1;Lc2;Lc3;Lc4;Lc5;Lc6;Lc7];
L = sparse(Lc(:,1),Lc(:,2),Lc(:,3),xnum*ynum,xnum*ynum);
Rc = Rc(:); Rc = cell2mat(Rc);
R = Rc(:,end);

% Solve matrix
S=L\R;

% Spar = S; Lpar = L; save('S_par','Spar'); save('L_par','Lpar');

% Reload solution
parfor i=1:ynum
    S = L\R;
    for j=1:xnum
        % Index for T
        itk=(j-1)*ynum+i;
        % Reload T
        tknew(i,j)=S(itk);
    end
end

% Computing residuals
rest = zeros(ynum,xnum); temp_rest1 = zeros(ynum,xnum); temp_rest2 = zeros(ynum,xnum);
parfor i=1:ynum
    temp_i = i;
    for j=1:xnum
        % Boundary conditions
        if (j==1 || j==xnum || temp_i==1 || temp_i==ynum)
            rest(i,j)=0;
        % Temperature equation RHO*CP*DT/Dt=-dqx/dx-dqy/dy+Ht
        else
            % Computing Current Temperature equation residual
            % Ht-DT/dt
            temp_rest1(i,j) = RT(temp_i,j)-rhocp(temp_i,j)*(tknew(temp_i,j)-tk(temp_i,j))/timestep;
            % -dqx/dx
            temp_rest2(i,j)=temp_rest1(i,j)+((kt(temp_i,j)+kt(temp_i,j+1))*(tknew(temp_i,j+1)-tknew(temp_i,j))/xstp(j)-(kt(temp_i,j-1)+kt(temp_i,j))*(tknew(temp_i,j)-tknew(temp_i,j-1))/xstp(j-1))/2/xstpc(j-1);
            % -dqy/dy
            rest(i,j)=temp_rest2(i,j)+((kt(temp_i,j)+kt(temp_i+1,j))*(tknew(temp_i+1,j)-tknew(temp_i,j))/ystp(temp_i)-(kt(temp_i-1,j)+kt(temp_i,j))*(tknew(temp_i,j)-tknew(temp_i-1,j))/ystp(temp_i-1))/2/ystpc(temp_i-1);
        end
    end
end
% restpar = rest;
% save('rest_par','restpar')
% a = DA_CANCELLARE;


