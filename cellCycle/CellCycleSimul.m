function varargout = CellCycleSimul(varargin)
% CELLCYCLESIMUL MATLAB code for CellCycleSimul.fig
%      CELLCYCLESIMUL, by itself, creates a new CELLCYCLESIMUL or raises the existing
%      singleton*.
%
%      H = CELLCYCLESIMUL returns the handle to a new CELLCYCLESIMUL or the handle to
%      the existing singleton*.
%
%      CELLCYCLESIMUL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLCYCLESIMUL.M with the given input arguments.
%
%      CELLCYCLESIMUL('Property','Value',...) creates a new CELLCYCLESIMUL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellCycleSimul_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellCycleSimul_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellCycleSimul

% Last Modified by GUIDE v2.5 29-Jul-2013 15:16:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellCycleSimul_OpeningFcn, ...
                   'gui_OutputFcn',  @CellCycleSimul_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before CellCycleSimul is made visible.
function CellCycleSimul_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellCycleSimul (see VARARGIN)

% Choose default command line output for CellCycleSimul
handles.output = hObject;
handles.liste = [];
handles.graphs=[];
handles.mode='';

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using CellCycleSimul.
if strcmp(get(hObject,'Visible'),'off')
    axes(handles.axes2);
    set(gca,'xtick',[],'ytick',[])
    cla;
end

% UIWAIT makes CellCycleSimul wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CellCycleSimul_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in submit.
function submit_Callback(hObject, eventdata, handles)
% hObject    handle to submit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a_axes=findall(0,'type','axes');
arrayfun(@cla,a_axes); %try not to change axe3 and axe4 which contains the image

 %mode= modes{get(handles.mode, 'Value')};
% divs=get(handles.enablediv, 'String');
% div= divs{get(handles.enablediv, 'Value')}

%% Get G1/S parameters
g1trans= str2double(get(handles.G1trans, 'String'));
g1trad= str2double(get(handles.G1trad, 'String'));
g1darn= str2double(get(handles.G1darn, 'String'));
g1dprot= str2double(get(handles.G1dprot, 'String'));
g1ig= str2double(get(handles.G1in, 'String'));
g1ag= str2double(get(handles.G1act, 'String'));
g1time= str2double(get(handles.G1time, 'String'));

%% Get G2 parameters
g2trans= str2double(get(handles.G2trans, 'String'));
g2trad= str2double(get(handles.G2trad, 'String'));
g2darn= str2double(get(handles.G2darn, 'String'));
g2dprot= str2double(get(handles.G2dprot, 'String'));
g2ig= str2double(get(handles.G2in, 'String'));
g2ag= str2double(get(handles.G2act, 'String'));
g2time= str2double(get(handles.G2time, 'String'));

%% Get M parameters
mtrans= str2double(get(handles.Mtrans, 'String'));
mtrad= str2double(get(handles.Mtrad, 'String'));
mdarn= str2double(get(handles.Mdarn, 'String'));
mdprot= str2double(get(handles.Mdprot, 'String'));
mig= str2double(get(handles.Min, 'String'));
mag= str2double(get(handles.Mact, 'String'));
mtime= str2double(get(handles.Mtime, 'String'));

%% Global parameters
arnnumber= str2double(get(handles.initARN, 'String'));
protnumber= str2double(get(handles.initPROT, 'String'));
time= str2double(get(handles.timer, 'String'));
ncycle= str2double(get(handles.ncycle, 'String'));

g1time=round(g1time*time/100);
g2time=round(g2time*time/100);
mtime=round(mtime*time/100);

g1 = [g1trans, g1trad, g1darn, g1dprot, g1ig, g1ag, g1time];
g2 = [g2trans, g2trad, g2darn, g2dprot, g2ig, g2ag, g2time];
m =  [mtrans, mtrad, mdarn, mdprot, mig, mag, mtime ];
if ((g1(5)~=0 || g1(6) ~=0) && (g2(5)~=0 || (g2(6)~=0)) && (m(5)~=0 || m(6)~=0))
    set(handles.mode, 'String', 'Burst');
else
    set(handles.mode, 'String','Cont');
end
[liste, handle]=CCSimulator(g1, g2, m, arnnumber, protnumber, time, ncycle);
% handles.liste= liste;
%  for i=1:numel(handle)
%      handles.graphs(i)=handle{i};
%  end
%  guidata(hObject, handles);
 
% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Are you sure you want to close this program?',...
    'Close Request','Yes','No','Yes');
switch selection
    case 'Yes',
        close all;
    case 'No'
        return
end %switch

% --- Executes on selection change in mode.
function mode_Callback(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mode


% --- Executes during object creation, after setting all properties.
function mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trans_Callback(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trans as text
%        str2double(get(hObject,'String')) returns contents of trans as a double


% --- Executes during object creation, after setting all properties.
function trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trad_Callback(hObject, eventdata, handles)
% hObject    handle to trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trad as text
%        str2double(get(hObject,'String')) returns contents of trad as a double


% --- Executes during object creation, after setting all properties.
function trad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function darn_Callback(hObject, eventdata, handles)
% hObject    handle to darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of darn as text
%        str2double(get(hObject,'String')) returns contents of darn as a double


% --- Executes during object creation, after setting all properties.
function darn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dprot_Callback(hObject, eventdata, handles)
% hObject    handle to dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dprot as text
%        str2double(get(hObject,'String')) returns contents of dprot as a double


% --- Executes during object creation, after setting all properties.
function dprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uiputfile(  {'untitle1.txt','TXT (*.txt)';'output.png','PNG (*.png)';'output.jpg','JPEG (*.jpg)';'output.tiff', 'Tiff (*.tiff)';'*.*', 'All file (*.*)'}, 'Save all as');
if(filterindex==1)
    file=fullfile(pathname, filename);
    dlmwrite(file, handles.liste.arn,'delimiter','\t');
    if isfield(handles.liste, 'prot');
        dlmwrite(strcat(file(1:end-4),'_p.txt') , handles.liste.prot,'delimiter','\t');
    end
elseif (filterindex~=0)
    imsave(filename, pathname, filterindex, handles);
end

% --------------------------------------------------------------------
function rotate_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%imrotate? 


% --- Executes on slider movement.
function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function timeslider_Callback(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.timer, 'String', num2str(round(get(hObject, 'Value'))));

% --- Executes during object creation, after setting all properties.
function timeslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function timer_Callback(hObject, eventdata, handles)
% hObject    handle to timer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timer as text
%        str2double(get(hObject,'String')) returns contents of timer as a double

min=get(handles.timeslider,'Min'); max=get(handles.timeslider,'Max');
val= str2double(get(hObject,'String'));
if(min<=val & val<=max)
    set(handles.timeslider, 'Value', val);
else
    warndlg('Warning, value exceeds bounds, it''ll be set to 10^3', 'Value');
    set(handles.timer, 'String', num2str(1000));
    set(handles.timeslider, 'Value', 1000);
end

% --- Executes during object creation, after setting all properties.
function timer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uiputfile(  {'output.png','PNG (*.png)';'output.jpg','JPEG (*.jpg)';'output.tiff', 'Tiff (*.tiff)';'*.*', 'All file (*.*)'}, 'Save image as');
imsave(filename, pathname, filterindex, handles);
%export_fig(h, fullfile(pathname, filename), '-png', '-m2.5');


function imsave(filename, pathname, filterindex, handles)
if exist('data.mat', 'file') && filterindex ~=0
    load data.mat;
    x=2; y=2;
    if exist('prot_array', 'var')
        x=3;
    end
    
    h=figure;
    set(h,'visible','off');
    subplot(x,y,[1 2]);
    %     plot(t_array, prom_array);
    %     xlabel('time (s)');
    %     ylabel('Promoter state.');
    %     set(gca, 'YTick', 0:1);
    %     set(gca,'YTickLabel',{'off', 'on'});
    
    t=1;
    for i=1:numel(prom_array)-1
        timearray(t)=t_array(i);
        timearray(t+1)=t_array(i);
        promoarray(t)=prom_array(i);
        if(prom_array(i)==prom_array(i+1))
            promoarray(t+1)=prom_array(i);
        else
            promoarray(t+1)=prom_array(i+1);
        end
        t=t+2;
    end
    plot(timearray, promoarray, '-r*');
    ylabel('Promoter state.');
    xlabel('time (s)');
    set(gca, 'YLim', [-1 2]);
    set(gca, 'YTick', [0 1]);
    set(gca,'YTickLabel',{'off', 'on'});
    
    
    subplot(x,y,3);
    plot(t_array,arn_array, 'g'); %mRNA time series
    xlabel('time (s)');
    ylabel('mRNA no.');
    
    subplot(x,y,4);
    mode=get(handles.mode, 'String');
    %
    if(isequal(mode , 'Burst'))
        hist(arn_array, max(arn_array)-min(arn_array)+1);
        set(gca, 'XLim', [min(arn_array) max(arn_array)]);
        set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
        legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)));
        
    else
        x= min(arn_array):max(arn_array);
        bar(x, hist(arn_array, max(arn_array)+1-min(arn_array))./ sum(hist(arn_array)), 'hist');%mRNA histogram
        set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
        legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)));
    end
    
    %%% Fitting a poisson distribution
    prob_distribution = fitdist(arn_array(:),'Poisson');
    hold on;
    probdf = pdf(prob_distribution,x);
    plot(x,probdf,'LineWidth',2, 'color', 'r');
    xlabel('mRNA no.');
    ylabel('percent');
    
    if(x==3)
        %%% Protein output
        subplot(x,y,5);
        plot(t_array,prot_array);
        xlabel('time (s)');
        ylabel('protein no.')
        
        subplot(x,y,6);
        hist(prot_array, max(prot_array)+1); %protein histogram
        xlabel('protein no.');
        ylabel('occurrence');
        legend(sprintf(' mean = %g \n stdev = %g',mean(prot_array), std(prot_array)))
    end
    
    saveas(h, fullfile(pathname, filename), 'png');
    close hidden;
end


% --------------------------------------------------------------------
function screenshot_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to screenshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uiputfile(  {'output.png','PNG (*.png)';'output.jpg','JPEG (*.jpg)';'output.tif', 'Tif (*.tif)';'output.pdf','PDF (*.pdf)';'*.*', 'All file (*.*)'}, 'Save image as');
switch(filterindex)
    case 1 
        export_fig(fullfile(pathname, filename), '-png', '-m2.5' );
    case 2 
        export_fig(fullfile(pathname, filename), '-jpg', '-m2.5');
    case 3 
        export_fig(fullfile(pathname, filename), '-tif', '-m2.5');
    case 4 
        export_fig(fullfile(pathname, filename), '-pdf', '-m2.5');
    case 5 
        export_fig(fullfile(pathname, filename));
end
 


% --------------------------------------------------------------------
function reset_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.G1trans, 'String', '0');
set(handles.G1trad, 'String', '0');
set(handles.G1darn, 'String', '0');
set(handles.G1dprot, 'String', '0');
set(handles.G1in, 'String', '0');
set(handles.G1act, 'String', '0');
set(handles.G1time, 'String', '0');

set(handles.G2trans, 'String', '0');
set(handles.G2trad, 'String', '0');
set(handles.G2darn, 'String', '0');
set(handles.G2dprot, 'String', '0');
set(handles.G2in, 'String', '0');
set(handles.G2act, 'String', '0');
set(handles.G2time, 'String', '0');

set(handles.Mtrans, 'String', '0');
set(handles.Mtrad, 'String', '0');
set(handles.Mdarn, 'String', '0');
set(handles.Mdprot, 'String', '0');
set(handles.Min, 'String', '0');
set(handles.Mact, 'String', '0');
set(handles.Mtime, 'String', '0');

set(handles.initPROT, 'String', '0');
set(handles.initARN, 'String', '0');



function celldiv_Callback(hObject, eventdata, handles)
% hObject    handle to celldiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of celldiv as text
%        str2double(get(hObject,'String')) returns contents of celldiv as a double


% --- Executes during object creation, after setting all properties.
function celldiv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to celldiv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in enablediv.
function enablediv_Callback(hObject, eventdata, handles)
% hObject    handle to enablediv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns enablediv contents as cell array
%        contents{get(hObject,'Value')} returns selected item from enablediv


% --- Executes during object creation, after setting all properties.
function enablediv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enablediv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncycle_Callback(hObject, eventdata, handles)
% hObject    handle to ncycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncycle as text
%        str2double(get(hObject,'String')) returns contents of ncycle as a double


% --- Executes during object creation, after setting all properties.
function ncycle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mtrans_Callback(hObject, eventdata, handles)
% hObject    handle to Mtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mtrans as text
%        str2double(get(hObject,'String')) returns contents of Mtrans as a double


% --- Executes during object creation, after setting all properties.
function Mtrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mdarn_Callback(hObject, eventdata, handles)
% hObject    handle to Mdarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mdarn as text
%        str2double(get(hObject,'String')) returns contents of Mdarn as a double


% --- Executes during object creation, after setting all properties.
function Mdarn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mdarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mdprot_Callback(hObject, eventdata, handles)
% hObject    handle to Mdprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mdprot as text
%        str2double(get(hObject,'String')) returns contents of Mdprot as a double


% --- Executes during object creation, after setting all properties.
function Mdprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mdprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Min_Callback(hObject, eventdata, handles)
% hObject    handle to Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Min as text
%        str2double(get(hObject,'String')) returns contents of Min as a double


% --- Executes during object creation, after setting all properties.
function Min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mact_Callback(hObject, eventdata, handles)
% hObject    handle to Mact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mact as text
%        str2double(get(hObject,'String')) returns contents of Mact as a double


% --- Executes during object creation, after setting all properties.
function Mact_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mtrad_Callback(hObject, eventdata, handles)
% hObject    handle to Mtrad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mtrad as text
%        str2double(get(hObject,'String')) returns contents of Mtrad as a double


% --- Executes during object creation, after setting all properties.
function Mtrad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mtrad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Minitarn_Callback(hObject, eventdata, handles)
% hObject    handle to Minitarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Minitarn as text
%        str2double(get(hObject,'String')) returns contents of Minitarn as a double


% --- Executes during object creation, after setting all properties.
function Minitarn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Minitarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Minitprot_Callback(hObject, eventdata, handles)
% hObject    handle to Minitprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Minitprot as text
%        str2double(get(hObject,'String')) returns contents of Minitprot as a double


% --- Executes during object creation, after setting all properties.
function Minitprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Minitprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mtime_Callback(hObject, eventdata, handles)
% hObject    handle to Mtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mtime as text
%        str2double(get(hObject,'String')) returns contents of Mtime as a double


% --- Executes during object creation, after setting all properties.
function Mtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2trans_Callback(hObject, eventdata, handles)
% hObject    handle to G2trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2trans as text
%        str2double(get(hObject,'String')) returns contents of G2trans as a double


% --- Executes during object creation, after setting all properties.
function G2trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2darn_Callback(hObject, eventdata, handles)
% hObject    handle to G2darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2darn as text
%        str2double(get(hObject,'String')) returns contents of G2darn as a double


% --- Executes during object creation, after setting all properties.
function G2darn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2dprot_Callback(hObject, eventdata, handles)
% hObject    handle to G2dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2dprot as text
%        str2double(get(hObject,'String')) returns contents of G2dprot as a double


% --- Executes during object creation, after setting all properties.
function G2dprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2in_Callback(hObject, eventdata, handles)
% hObject    handle to G2in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2in as text
%        str2double(get(hObject,'String')) returns contents of G2in as a double


% --- Executes during object creation, after setting all properties.
function G2in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2act_Callback(hObject, eventdata, handles)
% hObject    handle to G2act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2act as text
%        str2double(get(hObject,'String')) returns contents of G2act as a double


% --- Executes during object creation, after setting all properties.
function G2act_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2trad_Callback(hObject, eventdata, handles)
% hObject    handle to G2trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2trad as text
%        str2double(get(hObject,'String')) returns contents of G2trad as a double


% --- Executes during object creation, after setting all properties.
function G2trad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2initarn_Callback(hObject, eventdata, handles)
% hObject    handle to G2initarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2initarn as text
%        str2double(get(hObject,'String')) returns contents of G2initarn as a double


% --- Executes during object creation, after setting all properties.
function G2initarn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2initarn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2initprot_Callback(hObject, eventdata, handles)
% hObject    handle to G2initprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2initprot as text
%        str2double(get(hObject,'String')) returns contents of G2initprot as a double


% --- Executes during object creation, after setting all properties.
function G2initprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2initprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G2time_Callback(hObject, eventdata, handles)
% hObject    handle to G2time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G2time as text
%        str2double(get(hObject,'String')) returns contents of G2time as a double


% --- Executes during object creation, after setting all properties.
function G2time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G2time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1trans_Callback(hObject, eventdata, handles)
% hObject    handle to G1trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1trans as text
%        str2double(get(hObject,'String')) returns contents of G1trans as a double


% --- Executes during object creation, after setting all properties.
function G1trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1darn_Callback(hObject, eventdata, handles)
% hObject    handle to G1darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1darn as text
%        str2double(get(hObject,'String')) returns contents of G1darn as a double


% --- Executes during object creation, after setting all properties.
function G1darn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1dprot_Callback(hObject, eventdata, handles)
% hObject    handle to G1dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1dprot as text
%        str2double(get(hObject,'String')) returns contents of G1dprot as a double


% --- Executes during object creation, after setting all properties.
function G1dprot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1in_Callback(hObject, eventdata, handles)
% hObject    handle to G1in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1in as text
%        str2double(get(hObject,'String')) returns contents of G1in as a double


% --- Executes during object creation, after setting all properties.
function G1in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1act_Callback(hObject, eventdata, handles)
% hObject    handle to G1act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1act as text
%        str2double(get(hObject,'String')) returns contents of G1act as a double


% --- Executes during object creation, after setting all properties.
function G1act_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1act (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1trad_Callback(hObject, eventdata, handles)
% hObject    handle to G1trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1trad as text
%        str2double(get(hObject,'String')) returns contents of G1trad as a double


% --- Executes during object creation, after setting all properties.
function G1trad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initARN_Callback(hObject, eventdata, handles)
% hObject    handle to initARN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initARN as text
%        str2double(get(hObject,'String')) returns contents of initARN as a double


% --- Executes during object creation, after setting all properties.
function initARN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initARN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initPROT_Callback(hObject, eventdata, handles)
% hObject    handle to initPROT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initPROT as text
%        str2double(get(hObject,'String')) returns contents of initPROT as a double


% --- Executes during object creation, after setting all properties.
function initPROT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initPROT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function G1time_Callback(hObject, eventdata, handles)
% hObject    handle to G1time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of G1time as text
%        str2double(get(hObject,'String')) returns contents of G1time as a double


% --- Executes during object creation, after setting all properties.
function G1time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to G1time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
