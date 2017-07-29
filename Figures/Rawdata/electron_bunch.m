% Planned to show nice electron bunch with real dimensions in the paper
clear all;

%% Read in the image
image_folder = '\\cns\projects\HPLexp\Electrons\Experimental data\2016\LWFA\Lanex calibration\2016-05-10 Lanex calibration\DataI';
image_filename = 'Lanex_2016-05-11_00h-20m-00s_21018_original.png';
e_bunch = imread([image_folder '\' image_filename]);

% convert to double
e_bunch = double(e_bunch);

% normalize
e_bunch = e_bunch / max(e_bunch(:));

%Show raw image
figure
imagesc(e_bunch, 'Cdatamapping', 'scaled');

% Show ROI-image
pos_roi = [425 645 510 810];
e_bunch_roi = e_bunch(pos_roi(1):pos_roi(2), pos_roi(3):pos_roi(4));
figure
imagesc(e_bunch_roi, 'Cdatamapping', 'scaled');

%% Conversion into mm
conversion_x = 3.75*10^-3 * (362/53 -1);
conversion_y = 3.75*10^-3 * (362/53 -1);

%% Fitting & Plot
        
        % ---------User Input---------------------
        MdataSize_x = size(e_bunch_roi, 2); % Size of nxn data matrix
        MdataSize_y = size(e_bunch_roi, 1); % Size of nxn data matrix
        % parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
        x0 = [max(e_bunch_roi(:)),0,20,0,20,0]; %Inital guess parameters
        x = x0;
        noise = 0; % noise in % of centroid peak value (x(1))
        InterpolationMethod = 'linear'; % 'nearest','linear','spline','cubic'
        FitForOrientation = 0; % 0: fit for orientation. 1: do not fit for orientation

        % ---Generate centroid to be fitted--------------------------------------
        xin = x; 
        noise = noise/100 * x(1);
        X = linspace(-round(MdataSize_x/2)+1,round(MdataSize_x/2)-1,MdataSize_x);
        Y = linspace(-round(MdataSize_y/2)+1,round(MdataSize_y/2)-1,MdataSize_y);
        [X,Y] = meshgrid(X,Y); % generate high res grid for plot
        xdata = zeros(size(X,1),size(Y,2),2);
        xdata(:,:,1) = X;
        xdata(:,:,2) = Y;
        xdatahr = zeros(MdataSize_y,MdataSize_x,2);
        xdatahr(:,:,1) = X;
        xdatahr(:,:,2) = Y;
        %---Fit data
        Z = e_bunch_roi;

        % --- Fit---------------------
        if FitForOrientation == 0
            % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
            lb = [0,-MdataSize_y/2,0,-MdataSize_x/2,0,-pi/4];
            ub = [realmax('double'),MdataSize_y/2,(MdataSize_y/2)^2,MdataSize_x/2,(MdataSize_x/2)^2,pi/4];
            [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,Z,lb,ub);
        else
            x0 =x0(1:5);
            xin(6) = 0; 
            x =x(1:5);
            lb = [0,-MdataSize_y/2,0,-MdataSize_x/2,0];
            ub = [realmax('double'),MdataSize_y/2,(MdataSize_y/2)^2,MdataSize_x/2,(MdataSize_x/2)^2];
            [x,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub);
            x(6) = 0;
        end

        % -----Plot profiles----------------
        hf2 = figure;
        set(hf2, 'Position', [20 20 950 900])
        alpha(0)
        subplot('Position',[0.1 0.1 0.7 0.7]);
        imagesc(X(1,:)*conversion_x,Y(:,1)*conversion_y',Z);
        colormap('jet')
        xlabel('Length (mm)');
        ylabel('Length (mm)');
        %set(gca,'YDir','reverse');
        
        % -----Calculate cross sections-------------
        % generate points along horizontal axis
        m = -tan(x(6));% Point slope formula
        b = (-m*x(2) + x(4));
        xvh = -MdataSize_x/2:MdataSize_x/2;
        yvh = xvh*m + b;
        hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
        % generate points along vertical axis
        mrot = -m;
        brot = (mrot*x(4) - x(2));
        yvv = -MdataSize_y/2:MdataSize_y/2;
        xvv = yvv*mrot - brot;
        vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);

        hold on % Indicate major and minor axis on plot

        % % plot pints 
        % plot(xvh,yvh,'r.') 
        % plot(xvv,yvv,'g.')

        % plot lins 
        %plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
        %plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g') 
    
        hold off
        %axis([-MdataSize_x/2-0.5 MdataSize_x/2+0.5 -MdataSize_y/2-0.5 MdataSize_y/2+0.5])
        %

        ymin = - noise * x(1);
        ymax = x(1)*(1+noise);
        xdatafit = linspace(-MdataSize_x/2-0.5,MdataSize_x/2+0.5,300);
        hdatafit = x(1)*exp(-(xdatafit-x(2)).^2/(2*x(3)^2));
        vdatafit = x(1)*exp(-(xdatafit-x(4)).^2/(2*x(5)^2));
        subplot('Position',[0.1 0.8 0.7 0.12]);
        xposh = (xvh-x(2))/cos(x(6))+x(2);% correct for the longer diagonal if fi~=0
        plot(xposh,hPoints,'r','Linewidth',4);%'r.')%,xdatafit,hdatafit,'black')  
        
        axis([-MdataSize_x/2-0.5 MdataSize_x/2+0.5 ymin*1.2 ymax*1.2])
        set(gca,'XTick',zeros(1,0),'YTick',zeros(1,0));
        
        subplot('Position',[0.8 0.1 0.12 0.7]);
        xposv = (yvv-x(4))/cos(x(6))+x(4);% correct for the longer diagonal if fi~=0
        plot(vPoints,xposv,'g','Linewidth',4)%,vdatafit,xdatafit,'black')
        axis([ymin*1.1 ymax*1.1 -MdataSize_y/2-0.5 MdataSize_y/2+0.5])
        set(gca,'YDir','reverse','XTick',zeros(1,0),'YTick',zeros(1,0))
        figure(gcf) % bring current figure to front
        
        % Save the GaussFit
        cd('\\fwknas\kurz43\My Documents\Work\Projects\Lanex Calibration\Paper\Figures');
        saveas(gcf,'electron_bunch','epsc');  % Save the image