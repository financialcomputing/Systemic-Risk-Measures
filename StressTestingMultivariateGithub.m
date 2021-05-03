% code and data for the paper: 
% Stress testing and systemic risk measures using elliptical conditional multivariate probabilities
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3575512
% Tomaso Aste 03/05/2021
clear
close all

mod180 = @(x) x.*(x<=90)+(180-x).*(x>90);

load DataGithub.mat

O = cov(X0); %X0 are log-returns 

Angle = NaN(length(SectorNames));
DXY   = NaN(length(SectorNames));
TvXY  = NaN(length(SectorNames));
for s = 1:length(SectorNames)
    iX = find(strcmp(SectorNames(s),Sectors(:,3)));
    VaRX = quantile(X0(:,iX),1-q);
    OXX = O(iX,iX);
    for r = s:length(SectorNames)
        iY = find(strcmp(SectorNames(r),Sectors(:,3)));
        VaRY = quantile(X0(:,iY),1-q);
        OXY = O(iX,iY);
        OYX = O(iY,iX);
        OYY = O(iY,iY);
        %% columns | impact of 
        vLYX(r,s) = mean(OYX/OXX*(VaRX-mean(X0(:,iX)))'); %Y<-X s->r %the diagonal is the VaR (OYX/OXX=I)
        vLYX(s,r) = mean(OXY/OYY*(VaRY-mean(X0(:,iY)))'); %X<-Y r->s
        I(s,r) = 0.5*(sum(log2(eig(OXX)))+sum(log2(eig(OYY)))-sum(log2(eig(O))));
        %% angle 
        OYY_X = OYY - OYX/OXX*OXY ;
        if min(eig(OYY_X))<=0,fprintf('OYY_X not positive defined!\n')
        elseif ~isreal(eig(OYY_X)),fprintf('OYY_X not real !\n')
        else
        [u,e]=eig(OYY);
        [m,i]=max(max(e));
        uYY = u(:,i); 
        sYY = (m); 
        TvYY  = sum(log(diag(e)));
        [u,e]=eig(OYY_X);
        [m,i]=max(max(e));
        uYY_X = u(:,i); 
        sYY_X = (m); 
        TvYY_X  = sum(log(diag(e)));
        Angle(r,s) = acos(abs(uYY'*uYY_X))/pi*180; %Y <- X 
        DXY(r,s) = -(sYY_X-sYY)/sYY;
        TvXY(r,s)= (TvYY_X-TvYY);
        end
        OXX_Y = OXX - OXY/OYY*OYX ;
        if min(eig(OXX_Y))<=0,fprintf('OXX_Y not positive defined!\n')
        elseif ~isreal(eig(OXX_Y)),fprintf('OXX_Y not real !\n')
        else
        [u,e]=eig(OXX);
        [m,i]=max(max(e));
        uXX = u(:,i); 
        sXX = (m); 
        TvXX  = sum(log(diag(e)));
        [u,e]=eig(OXX_Y);
        [m,i]=max(max(e));
        uXX_Y = u(:,i); 
        sXX_Y = (m); 
        TvXX_Y  = sum(log(diag(e)));
        Angle(s,r) = acos(abs(uXX'*uXX_Y))/pi*180; %X <- Y  
        DXY(s,r) = -(sXX_Y-sXX)/sXX;
        TvXY(s,r)= (TvXX_Y-TvXX);        
        end
    end
    %% Mahalanobis factor
    dX(s) = (quantile(X0(:,iX),1-q)-mean(X0(:,iX)))/OXX*(quantile(X0(:,iX),1-q)-mean(X0(:,iX)))';
    pX(s) = length(iX);
    kY = setdiff([1:N],iX)';
    OYY = O(kY,kY);
    OXY = O(iX,kY);
    OXX = O(iX,iX);
    OYX = O(kY,iX);
    %% Impacted and impacting on/from the rest of the system
    vImpacted(s)  =  mean(OXY/OYY*(quantile(X0(:,kY),1-q)-mean(X0(:,kY)))'); % impact on s from all others
    vImpacting(s) =  mean(OYX/OXX*(quantile(X0(:,iX),1-q)-mean(X0(:,iX)))'); % impact of s on all others
end
%% Average losses 
% reading: 
% vLYX(r,c)  Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!! therefore y-axis -> x-axis in the figure 
figure
M=-(vLYX-diag(diag(nan(10))))';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
colormap('hot')
colorbar
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Imapct map')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)

%%
figure
plot([0 2],[0 2],'m')
hold on
plot(-vImpacted,-vImpacting,'oc','markersize',18) % impact of the other on a given sector
text(-vImpacted,-vImpacting,SectorNames,'fontname','times','fontsize',14)
xlim([0.02 .25])
ylim([0.04 .105])
set(gca,'fontsize',24)
xlabel('Average losses on the sector (VaR$_{0.95}$ stress)','interpreter','latex','fontsize',20)
ylabel('Average losses from the sector (Var$_{0.95}$ stress)','interpreter','latex','fontsize',20)

%%
figure
im=imagesc((I+I'-diag(diag(nan(10)))))
set(im,'AlphaData',~isnan((I+I'-diag(diag(nan(10))))))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan((I+I'-diag(diag(nan(10))))));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
%
title('Mutual Information')
set(gca,'fontsize',14,'fontname','times')


%% angle rotation  
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = mod180(Angle-diag(diag(Angle)))';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Principal axis rotation angle')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)

%% eigenvalue change , 
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = DXY';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Max eigenvalue change')
set(gca,'fontsize',14,'fontname','times')
ylabel('stress','interpreter','latex','fontsize',20)
xlabel('response','interpreter','latex','fontsize',20)


%% Total Variance
% reading: 
% Impact of c -> on r :: columns -> rows
% BUT THE FIGURE IS TRANSPOSED!!!!!
figure
M = TvXY';
im=imagesc(M)
set(im,'AlphaData',~isnan(M))
colormap('hot')
colorbar
hold on
[i,j]=find(isnan(M));
plot(i,j,'xb','markersize',30)
xticklabels(cellstr(SectorNames))
xtickangle(60)
yticklabels(cellstr(SectorNames))
title('Total Variance log-change ')
set(gca,'fontsize',14,'fontname','times')

%% Mahlanobis impact factor
figure
plot(dX./pX,'o-')
xticklabels(cellstr(SectorNames))
xtickangle(60)
title('Mahlanobis impact factor')
set(gca,'fontsize',14,'fontname','times')
xlabel('conditioning variables','interpreter','latex','fontsize',20)
ylabel('$d_{\mathbf x_q}/p_{\mathbf X}$','interpreter','latex','fontsize',20)
