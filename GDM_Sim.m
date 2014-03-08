function varargout = GDM_Sim(varargin)
% GDM_SIM MATLAB code for GDM_Sim.fig
%      GDM_SIM, by itself, creates a new GDM_SIM or raises the existing
%      singleton*.
%
%      H = GDM_SIM returns the handle to a new GDM_SIM or the handle to
%      the existing singleton*.
%
%      GDM_SIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GDM_SIM.M with the given input arguments.
%
%      GDM_SIM('Property','Value',...) creates a new GDM_SIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GDM_Sim_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GDM_Sim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GDM_Sim

% Last Modified by GUIDE v2.5 18-Jul-2013 11:42:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GDM_Sim_OpeningFcn, ...
                   'gui_OutputFcn',  @GDM_Sim_OutputFcn, ...
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

% --- Executes just before GDM_Sim is made visible.
function GDM_Sim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GDM_Sim (see VARARGIN)

% Choose default command line output for GDM_Sim
handles.output = hObject;
handles.liste = [];
handles.graphs=[];

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using GDM_Sim.
if strcmp(get(hObject,'Visible'),'off')
    axes(handles.axe3);
    %set(gca, 'visible', 'off');
    imshow(imread('burst.png'));
    axes(handles.axe4);
    %set(gca, 'visible', 'off');
    imshow(imread('continous.png'));
    axes(handles.axes2);
    set(gca,'xtick',[],'ytick',[])
    cla;
end

% UIWAIT makes GDM_Sim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GDM_Sim_OutputFcn(hObject, eventdata, handles)
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
arrayfun(@cla,a_axes(1:end-2)); %try not to change axe3 and axe4 which contains the image
modes=get(handles.mode, 'String');
mode= modes{get(handles.mode, 'Value')};
divs=get(handles.enablediv, 'String');
div= divs{get(handles.enablediv, 'Value')}
trans= str2double(get(handles.trans, 'String'));
trad= str2double(get(handles.trad, 'String'));
darn= str2double(get(handles.darn, 'String'));
dprot= str2double(get(handles.dprot, 'String'));
ig= str2double(get(handles.ig, 'String'));
ag= str2double(get(handles.ag, 'String'));
arnnumber= str2double(get(handles.arnnumber, 'String'));
protnumber= str2double(get(handles.protnumber, 'String'));
time= str2double(get(handles.timer, 'String'));
celldiv= str2double(get(handles.celldiv, 'String'));
[liste, handle]=GSimulator(mode, time, trans, darn,arnnumber, trad, dprot, protnumber, ag, ig, div, celldiv);
handles.liste= liste;
 for i=1:numel(handle)
     handles.graphs(i)=handle{i};
 end
 guidata(hObject, handles);
 
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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
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



function ip_Callback(hObject, eventdata, handles)
% hObject    handle to ip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ip as text
%        str2double(get(hObject,'String')) returns contents of ip as a double


% --- Executes during object creation, after setting all properties.
function ip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ap_Callback(hObject, eventdata, handles)
% hObject    handle to ap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ap as text
%        str2double(get(hObject,'String')) returns contents of ap as a double


% --- Executes during object creation, after setting all properties.
function ap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of darn as text
%        str2double(get(hObject,'String')) returns contents of darn as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dprot as text
%        str2double(get(hObject,'String')) returns contents of dprot as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dprot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ig_Callback(hObject, eventdata, handles)
% hObject    handle to ig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ig as text
%        str2double(get(hObject,'String')) returns contents of ig as a double


% --- Executes during object creation, after setting all properties.
function ig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ag_Callback(hObject, eventdata, handles)
% hObject    handle to ag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ag as text
%        str2double(get(hObject,'String')) returns contents of ag as a double


% --- Executes during object creation, after setting all properties.
function ag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trad as text
%        str2double(get(hObject,'String')) returns contents of trad as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function arnnumber_Callback(hObject, eventdata, handles)
% hObject    handle to arnnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of arnnumber as text
%        str2double(get(hObject,'String')) returns contents of arnnumber as a double


% --- Executes during object creation, after setting all properties.
function arnnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to arnnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function protnumber_Callback(hObject, eventdata, handles)
% hObject    handle to protnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of protnumber as text
%        str2double(get(hObject,'String')) returns contents of protnumber as a double


% --- Executes during object creation, after setting all properties.
function protnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to protnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


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


% --- Executes during object creation, after setting all properties.
function axe3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axe3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axe3


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
    modes=get(handles.mode, 'String');
    mode= modes{get(handles.mode, 'Value')};
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
set(handles.trans, 'String', '0');
set(handles.trad, 'String', '0');
set(handles.darn, 'String', '0');
set(handles.dprot, 'String', '0');
set(handles.ig, 'String', '0');
set(handles.ag, 'String', '0');
set(handles.arnnumber, 'String', '0');
set(handles.protnumber, 'String', '0');
set(handles.celldiv, 'String', '0');



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
