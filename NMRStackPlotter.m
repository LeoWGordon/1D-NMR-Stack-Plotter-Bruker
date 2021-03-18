% Choose the experiment FOLDER of the experiment and fill the proc. no box.
% GUI to plot 1D NMR spectra stacked
% After the top spectrum is plotted, each subsequent one will plot below
% and the axes will update to fit
% 1H, 13C, 27Al, and 14N are available with the current limits, but easy to
% change.
%%%%%%%
%% Version History (Latest to oldest)
%% Version 2.1 (Dec 29th 2020)
% Version Changes:-
% - Added a method to plot a zoomed portion of a line through mouse clicks
% - Adds a label to the zoomed line to state level of magnification
%%%%%%%
%% Version 2.0 (Nov 24th 2020)
% A whole version update due to the major changes under the hood
% Version changes:-
% - Changed the 1D import function completely, reads the spectrum file
% directly rather than needing the convbin2asc au program.
% - Added a file type selection for saving (pdf, png, fig)
% - Added manual axis limits (which update the defaults
%%%%%%%
%% Version 1.1
% Version changes:
% -Added radio button to choose colours or all black
% -Added option to remove the most recent spectrum
%%%%%%%

%%
function NMRStackPlotter
% Create figure
h.f(1) = figure('units','normalized','position',[0.05,0.15,0.3,0.6],...
             'toolbar','none','menu','none');
h.c(1) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.25,0.9,0.55,0.05],'String','Choose Experiment Folder','FontSize',14); %1st file input
h.c(2) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.8,0.76,0.15,0.05],'String','Nucleus','FontSize',14); %F1 nucleus string
h.c(3) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.8,0.73,0.15,0.05],'String',{'<HTML><sup>1</sup>H','<HTML><sup>13</sup>C','<HTML><sup>27</sup>Al','<HTML><sup>14</sup>N'},'Callback',@updatelimits); %F2 dropdown
h.c(4) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.25,0.65,0.55,0.05],'String','Input Subsequent Experiments','FontSize',14); %1st file input
h.c(5) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.8,0.935,0.145,0.05],'String','Proc. No','FontSize',14); % F1 nucleus string
h.c(6) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.8,0.9,0.145,0.05],'String','1'); % Proc. no
h.c(7) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.8,0.685,0.145,0.05],'String','Proc. No','FontSize',14); % F1 nucleus string
h.c(8) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.8,0.65,0.145,0.05],'String','1'); % Proc. no
h.c(9) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.15,0.735,0.25,0.035],'String','','FontSize',14); % Lower limit Value
h.c(10) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.4,0.735,0.25,0.035],'String','','FontSize',14); % Upper limit Value
h.c(11) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.15,0.775,0.25,0.035],'String','Lower Limit','FontSize',14); % Lower limit Text
h.c(12) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.4,0.775,0.25,0.035],'String','Upper Limit','FontSize',14); % Upper limit Text

% magnify range of a chosen line (to do after all lines are plotted?)
            
% populate h.c(9) and h.c(10) based on nucleus choice
        nucvalintro = get(h.c(3), 'Value');

            if nucvalintro == 1
                lowerlim = -2;
                upperlim = 18;
            elseif nucvalintro == 2
                lowerlim = 0;
                upperlim = 200;
            elseif nucvalintro == 3
                lowerlim = -60;
                upperlim = 140;
            elseif nucvalintro == 4
                lowerlim = -200;
                upperlim = 200;
            end
        set(h.c(9),'String',num2str(lowerlim));
        set(h.c(10),'String',num2str(upperlim))
        
    function updatelimits(varargin)
        nucvalupdate = get(h.c(3), 'Value');

            if nucvalupdate == 1
                lowerlimupdate = -2;
                upperlimupdate = 18;
            elseif nucvalupdate == 2
                lowerlimupdate = 0;
                upperlimupdate = 200;
            elseif nucvalupdate == 3
                lowerlimupdate = -60;
                upperlimupdate = 140;
            elseif nucvalupdate == 4
                lowerlimupdate = -200;
                upperlimupdate = 200;
            end
        set(h.c(9),'String',num2str(lowerlimupdate))
        set(h.c(10),'String',num2str(upperlimupdate))
    end
            
% add in some radio buttons for choosing colours in the subsequent plot
% section
bg = uibuttongroup('Visible','off','units','normalized',...
                'position',[0.1,0.435,0.8,0.175],'BorderType','none');
h.r(1) = uicontrol(bg,'Style','radiobutton','units','normalized',...
                'position',[0.2,0.425,0.3,0.175],'string','Colours','FontSize',14);            
h.r(2) = uicontrol(bg,'Style','radiobutton','units','normalized',...
                'position',[0.6,0.425,0.3,0.175],'string','Monochrome','FontSize',14);

bg.Visible = 'on';

h.c(13) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.605,0.445,0.2,0.05],'String','File Type','FontSize',14); % File type to save
h.c(14) = uicontrol('Style','edit','Units','normalized',...
                'Position',[0.05,0.425,0.55,0.05],'String','Input Save File Name','FontSize',14); %Save name input
h.c(15) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.605,0.41,0.2,0.05],'String',{'pdf','png','fig'}); % File type to save
% Magnification Controls
h.c(16) = uicontrol('Style','text','Units','normalized',...
                'Position',[0.55,0.165,0.2,0.035],'String','Magnification Level','FontSize',14); % Magnification Box Text
h.c(17) = uicontrol('Style','popupmenu','Units','normalized',...
                'Position',[0.75,0.15,0.15,0.05],'String',{'2x','4x','8x','16x','32x'},'Callback',@magnumber); % Magnification dropdown

            magnumber_initial = 2^h.c(17).Value;
            
            DATA.magnumber = magnumber_initial;
            
% Create Pushbuttons
h.p(1) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.825,0.9,0.075],'string','Initial Plot','FontSize',14,...
                'callback',@p_load);
h.p(2) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.575,0.9,0.075],'string','Plot','FontSize',14,...
                'callback',@p_plot);
h.p(3) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.8,0.425,0.15,0.05],'string','Save','FontSize',14,...
                'callback',@p_save);
h.p(4) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.275,0.9,0.05],'string','Clear Figures','FontSize',14,...
                'callback',@p_clear);
h.p(5) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.9,0.2,0.05],'string','Browse Data','FontSize',14,...
                'callback',@p_browse_initial); % initial plot browse
h.p(6) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.65,0.2,0.05],'string','Browse Data','FontSize',14,...
                'callback',@p_browse_subsequent); % subsequent plot browse
h.p(7) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.05,0.35,0.9,0.05],'string','Remove Last Plot','FontSize',14,...
                'callback',@p_clear_last);
h.p(8) = uicontrol('style','pushbutton','units','normalized',...
                'position',[0.5,0.085,0.45,0.05],'string','Plot Magnified Section','FontSize',14,...
                'callback',@p_mag);            
% This creates the 'background' axes
ha = axes('units','normalized', ...
            'position',[0.05 0.03 0.4 0.22]);
% Move the background axes to the bottom
uistack(ha,'bottom');
% Load in a background image and display it using the correct colors
I = imread('./Messinger_Logo_Design.png');
hi = imagesc(I);
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
set(ha,'handlevisibility','off', ...
            'visible','off')  
    
% Pushbutton callbacks
%%%%%%        
    function p_browse_initial(varargin)
        
        f2name = uigetdir('~/','Choose the folder of the experiment');
        if f2name == 0
            set(h.c(1),'String','Input Experiment Folder')
        else
            set(h.c(1),'String',f2name)
        end

    end
%%%%%%

%%%%%%
    function p_browse_subsequent(varargin)
        
        f1name = uigetdir('~','Choose the folder of the experiment');
        if f1name == 0
            set(h.c(4),'String','Input Experiment Folder')
        else
            set(h.c(4),'String',f1name)
        end

    end
%%%%%%
       
%%%%%%
    function p_load(varargin)
            % Load first spectrum, creates first axis
            nucval = get(h.c(3), 'Value');
            
            filename = h.c(1).String;
            procinitial = str2double(h.c(6).String);
            [x1,y1,~]=brukimp1d(filename,procinitial); % change this for new import fn
            
            h.f(2) = figure;
            H = plot(x1,y1,'.-k','LineWidth',2);
            
            if nucval == 1
                xticks(-200:2:300)
                xlabel('^{1}H Chemical Shift (ppm)')
            elseif nucval == 2
                xticks(-200:20:300)
                xlabel('^{13}C Chemical Shift (ppm)')
            elseif nucval == 3
                xticks(-200:20:300)
                xlabel('^{27}Al Shift (ppm)')
            elseif nucval == 4
                xticks(-200:20:300)
                xlabel('^{14}N Shift (ppm)')
            end
            
            xlower = str2double(h.c(9).String);
            xupper = str2double(h.c(10).String);
            axnum = [xlower xupper 0 1];
            axis(axnum);            
            ax = gca;
            ax.FontName = 'Helvetica';
            ax.FontWeight = 'bold';
            ax.LineWidth = 2;
            legend off
            ax.FontSize = 20;
            ax.XDir = 'reverse';
            ax.TickDir = 'out';
            ax.YAxis.Visible = 'off';
            ax.XGrid = 'off';
            ax.XMinorGrid = 'off';
            ax.XMinorTick = 'on';
            ax.YGrid = 'off';
            axis square
            box off
            ax.XDir = 'reverse';
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperUnits', 'inches');
            set(gcf, 'PaperPosition', [0 0 10 11]);
            
            hold on
                    
            DATA.ax = ax;
            DATA.axnum = axnum;
            DATA.H1 = H;

            guidata(h.f(1),DATA);
    end
%%%%%%

%%%%%%
    function p_plot(varargin)
            
        DATA = guidata(h.f(1));
        ax = DATA.ax;
        axnum = DATA.axnum;
        procafter = str2double(h.c(6).String);
        
        [xn,yn,~] = brukimp1d(h.c(4).String, procafter);
        n = length(findall(h.f(2),'type','line'));
        i=n;
        
        if h.r(1).Value == 1

            cmapfig = colormap([1 0 0; 0 0 1; 1 0 1; 0 1 0; 0 1 1; 1 1 0; 0 0 0]);
            pc = cmapfig(n,:);

            H(i) = plot(xn,yn-n,'.-','Color',pc,'LineWidth',2,'Parent',ax);
        elseif h.r(2).Value == 1
            H = plot(xn,yn-n,'.-k','LineWidth',2,'Parent',ax);
        end

        axnum(3) = -n;
        
        ax.YLim = axnum(3:4);
        
        set(h.c(4),'String','')
        
        DATA.H(i) = H(i);

        guidata(h.f(1),DATA)
    end
%%%%%%

%%%%%%
    function p_save(varargin)
        savenmrfig(h.c(14).String)
                function savenmrfig(name)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Here you can make the savelocation variable to suit 
                    % your computer. Change to the path of the desired folder
                    savelocation = "~";
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % maybe need to assign numbers because dropdown doesnt
                    % call the strings or something, 1=pdf,2=png,3=fig.
                    type = h.c(15).Value;
                    
                    orient(h.f(2),'landscape')

                    if type==1
                        filenamepdf=strcat(string(name), '.pdf');
%                         orient(h.f(2),'landscape')
                        h.f(2).Renderer='painters'; % ensures the figure is actually saved as a vector graphic
                        print(h.f(2),filenamepdf,'-dpdf','-fillpage')
                        movefile(filenamepdf,savelocation)
                    elseif type==2
                        filenamepng=strcat(string(name), '.png');
                        print(h.f(2),filenamepng,'-dpng','-r1000')
                        movefile(filenamepng,savelocation)
                    elseif type==3
                        figname=strcat(string(name), '.fig');
                        saveas(h.f(2),name)
                        movefile(figname,savelocation)
                    end
                    
                end
    end
%%%%%%

%%%%%%
    function p_clear(varargin)
        set(h.f(1), 'HandleVisibility', 'off')
        close all;
        set(h.f(1), 'HandleVisibility', 'on')
        set(h.c(14),'String','Input Save File Name')
        set(h.c(1),'String','')
    end
%%%%%%

%%%%%%
    function p_clear_last(varargin)
        % Import figure properties
        DATA = guidata(h.f(1));
        ax = DATA.ax;
        axnum = DATA.axnum;
        
        % Find number of lines and generate handle for most recent
        n = length(findall(h.f(2),'type','line'));
        l = findobj('type','line');
        % Remove last line
        delete(l(1))
        % Fix axes
        axnum(3) = -(n-2);
        ax.YLim = axnum(3:4);
    end
%%%%%%
    function magnumber(varargin)
            magupdate = h.c(17).Value;

            magnumber = 2^magupdate;
            
            DATA.magnumber = magnumber;
            guidata(h.f(1),DATA);
    end

    function DATA = p_mag(varargin)
        
        DATA = guidata(h.f(1));
        H1 = DATA.H1;
%         magnification = DATA.magnumber;
        try 
            H =  DATA.H;  
            set(H, 'ButtonDownFcn', {@LineSelected, H})
        catch
        end     
       
        set(H1, 'ButtonDownFcn', {@LineSelected, H})
        
        
        function LineSelected(ObjectH, EventData, H)
            
            DATA = guidata(h.f(1));
            ax = DATA.ax;
            magnification = DATA.magnumber;

            [magx, ~] = ginput_ax(2);
            xdata = get(ObjectH, 'XData')';
            ydata = get(ObjectH, 'YData')';
            colour = get(ObjectH, 'Color');
            
            if magx(2)>magx(1)
                istart = find(xdata>=magx(2),1,'last'); %gets higher ppm value index
                iend = find(xdata>=magx(1),1,'last'); % gets lower ppm value index
            else
                istart = find(xdata>=magx(1),1,'last'); %gets higher ppm value index
                iend = find(xdata>=magx(2),1,'last'); % gets lower ppm value index
            end
            
            baseline = ydata(istart);
            
            xmag = xdata(istart:iend);
            ymag = ydata(istart:iend).*magnification;
            
            magbaseline = min(ymag);
            ymag = (ymag)+baseline-magbaseline+0.1;
            
            plot(ax,xmag,ymag,'.-','Color',colour)
            
            % adds magnification label
            text(xmag(end)-0.5,ymag(end),strcat(['\times',num2str(magnification)]),'FontSize',12,'FontWeight','bold')
            
        end

    end

% Functions called:
%%%%%%
    function [xppm, Spectrum, Params] = brukimp1d(fileinput, proc_no)

        input_file = strcat(fileinput,'/pdata/',string(proc_no),'/1r');
        Acqus_1D_H_File_Path = string([fileinput '/acqus']);
        Procs_1D_H_File_Path = strcat(fileinput, '/pdata/', string(proc_no), '/procs');

        fid = fopen(input_file, 'rb');
        if fid < 1
            error('File not found %s\n', input_file);
        else
            Spectrum_1D = fread(fid, 'int');
        end
        fclose(fid);

        fid_aqus = fopen(Acqus_1D_H_File_Path, 'r');
        fid_procs = fopen(Procs_1D_H_File_Path, 'r');
        if fid_aqus < 1 || fid_procs < 1
            fclose(fid_aqus);
            fclose(fid_procs);
            error('Could not open %s or %s\n', Acqus_1D_H_File_Path, Procs_1D_H_File_Path);
        else
            [H_OBS, H_CAR, H_Error_aqus] = Get_Bruker_Info_1D_Acqus(fid_aqus);
            [H_SF, H_SW, H_Length, H_Error_proc] = Get_Bruker_Info_1D_Procs(fid_procs);
            if ~isempty(H_Error_aqus) || ~isempty(H_Error_proc)
                fclose(fid_aqus);
                fclose(fid_procs);
                error('Something went wrong with the params in %s or %s\n', Acqus_1D_H_File_Path, Procs_1D_H_File_Path);
            end
        end
        fclose(fid_aqus);
        fclose(fid_procs);

        Params.xOBS = H_OBS;
        Params.xCAR = H_CAR;
        Params.xSW = H_SW;
        Params.xSF = H_SF;

        Spectrum = Spectrum_1D;
        Spectrum = Spectrum-min(Spectrum);
        Spectrum = Spectrum./max(Spectrum);

        xaxcen=Params.xCAR*Params.xOBS-((Params.xSF-Params.xOBS)*1000000); % why don't I just read in SFO1? not same?
        xaxmin=xaxcen-Params.xSW/2;
        xaxmax=xaxcen+Params.xSW/2;
        xaxlen=(xaxmax-xaxmin)/(H_Length-1);
        xaxhz=(xaxmin:xaxlen:xaxmax);
        xppm=xaxhz/Params.xOBS;
        xppm=sort(xppm,'descend');
        
%%%%%%

%%%%%%
        function [SF, SW, Length, Error] = Get_Bruker_Info_1D_Procs(fid)

            SW = 0;
            Length = 0;
            SF = 0;

            tline = fgetl(fid);
            Satisfied = false;
            while ~Satisfied
                    if ~isempty(strfind(tline, '##$SW_p= '))
                        tline = strrep(tline, '##$SW_p= ', '');
                        SW = str2double(tline);
                    end
                    if ~isempty(strfind(tline, '##$SI= '))
                        tline = strrep(tline, '##$SI= ', '');        
                        Length = str2double(tline);
                    end
                    if ~isempty(strfind(tline, '##$SF= '))
                        tline = strrep(tline, '##$SF= ', '');
                        SF = str2double(tline);
                    end
                tline = fgetl(fid);
                    if ~ischar(tline) || (SW~=0 && Length~= 0 && SF~=0)
                        Satisfied = true;
                    end
            end


                if (SW~=0 && Length ~= 0)
                    Error = '';
                else
                    Error = 'Could not find all the parameters from the aqcus file';
                end
        end
%%%%%%

%%%%%%
        function [OBS, CAR, Error] = Get_Bruker_Info_1D_Acqus(fid)

            OBS = nan;
            CAR = nan;
            O1 = nan;

            tline = fgetl(fid);
            Satisfied = false;
            while ~Satisfied
                if ~isempty(strfind(tline, '##$O1= '))
                    tline = strrep(tline, '##$O1= ', '');
                    O1 = str2double(tline);
                end
                if ~isempty(strfind(tline, '##$BF1= '))
                    tline = strrep(tline, '##$BF1= ', '');
                    OBS = str2double(tline);
                end
                tline = fgetl(fid);
                if ~isnan(OBS) && ~isnan(O1)
                    Satisfied = true;
                end
            end

            if (OBS~=0 && O1~=0)
                CAR = O1/OBS;
                Error = '';
            elseif O1==0
                CAR = 1e-30;
                Error = '';
            else
                Error = 'Could not find all the parameters from the aqcus file';
            end
        end
    end
%%%%%%

%%%%%%
    function varargout = ginput_ax(n)

    k = 0;
    xy = zeros(n,2);
    hf = h.f(2);
    ha = hf.CurrentAxes;
%     hf = get(ha,'parent');
    figure(hf);
    set(hf,'WindowButtonMotionFcn',@changepointer)
    set(ha,'ButtonDownFcn',@getpoints)
    hp = get(ha,'children');
    ht = get(hp,'hittest');
    set(hp,'hittest','off')
    axlim = get(ha,'Position');
    fglim = get(hf,'Position');
    x1 = axlim(1)*fglim(3) + fglim(1);
    x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
    y1 = axlim(2)*fglim(4) + fglim(2);
    y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
    waitfor(hf,'WindowButtonMotionFcn',[])
    if iscell(ht)
        for jj=1:length(ht)
            set(hp(jj),'hittest',ht{jj})
        end
    else
        set(hp,'hittest',ht)
    end
    if nargout==2
        varargout{1} = xy(:,1);
        varargout{2} = xy(:,2);
    else
        varargout{1} = xy;
    end
          function changepointer(~,~)
              pntr = get(0,'PointerLocation');
              if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
                  set(hf,'Pointer','crosshair')
              else
                  set(hf,'Pointer','arrow')
              end
          end
          function getpoints(hObj,~,~)
              cp = get(hObj,'CurrentPoint');
              k = k+1;
              xy(k,:) = cp(1,1:2);
              if k==n
                  set(hf,'Pointer','arrow')
                  set(hf,'WindowButtonMotionFcn',[])
                  set(ha,'ButtonDownFcn',[])
              end
          end
    end
%%%%%%
end