function varargout = behavior_correct_speed(varargin)
% BEHAVIOR_CORRECT_SPEED MATLAB code for behavior_correct_speed.fig
%      BEHAVIOR_CORRECT_SPEED, by itself, creates a new BEHAVIOR_CORRECT_SPEED or raises the existing
%      singleton*.
%
%      H = BEHAVIOR_CORRECT_SPEED returns the handle to a new BEHAVIOR_CORRECT_SPEED or the handle to
%      the existing singleton*.
%
%      BEHAVIOR_CORRECT_SPEED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEHAVIOR_CORRECT_SPEED.M with the given input arguments.
%
%      BEHAVIOR_CORRECT_SPEED('Property','Value',...) creates a new BEHAVIOR_CORRECT_SPEED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before behavior_correct_speed_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to behavior_correct_speed_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help behavior_correct_speed

% Last Modified by GUIDE v2.5 16-Aug-2019 13:43:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @behavior_correct_speed_OpeningFcn, ...
                   'gui_OutputFcn',  @behavior_correct_speed_OutputFcn, ...
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

% --- Executes just before behavior_correct_speed is made visible.
function behavior_correct_speed_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to behavior_correct_speed (see VARARGIN)

% Choose default command line output for behavior_correct_speed
handles.output = hObject;

handles.img_stack=varargin{1};
handles.filename=varargin{2};
handles.istart=varargin{3};
handles.iend=varargin{4};
handles.dorsal=varargin{5};
handles.ventral=varargin{6};
handles.centerline=varargin{7};
handles.velant=varargin{8};
handles.velpos=varargin{9};
handles.smallarea=varargin{10};

if handles.smallarea==1
    handles.img_stack(:,1) = cellfun(@(x) x(end/4:end/4*3, end/4:end/4*3), handles.img_stack(:,1), 'UniformOutput', false);
end

[height,width]=size(handles.img_stack{1,1});
handles.image_width=width;
handles.image_height=height;
handles.image_depth=handles.iend-handles.istart+1;
sections=length(handles.img_stack);
handles.sections=sections;
handles.framenumber = 1;

axes(handles.Velocity);
hold off;
% plot(handles.velant,'g');
% hold on;
plot(handles.velpos,'k'); hold on;
plot([1 handles.image_depth], [0 0], 'r');
handles.med = median(abs([handles.velant handles.velpos]), 'all');
handles.mag = 5;
ylim([-handles.mag*handles.med handles.mag*handles.med]);
plot([handles.framenumber handles.framenumber], ylim, 'm');

axes(handles.Recording); hold off;
img=imagesc(handles.img_stack{1,1});
colormap(gray);
hold on;
i=handles.istart;
plot(handles.dorsal{2*i-1,1}, handles.dorsal{2*i,1}, 'r');
plot(handles.ventral{2*i-1,1}, handles.ventral{2*i,1}, 'r');
plot(handles.centerline(:, 2*i-1), handles.centerline(:, 2*i), 'w', 'linewidth', 1.5);
plot(handles.centerline(1,2*i-1), handles.centerline(1, 2*i), 'og', 'markersize', 4, 'markerfacecolor', 'g');
plot(handles.centerline(end/2,2*i-1), handles.centerline(end/2, 2*i), 'ok', 'markersize', 4, 'markerfacecolor', 'k');
plot(handles.centerline(end,2*i-1), handles.centerline(end, 2*i), 'oy', 'markersize', 4, 'markerfacecolor', 'y');
text(10,10,['Frame #' num2str(handles.framenumber)],'color', 'k');
set(handles.Recording, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0], 'Units', 'pixels', 'DataAspectRatio', [1 1 1]);
set(img,'ButtonDownFcn', 'proof_reading(''ButtonDown_Callback'',gcbo,[],guidata(gcbo))');

min_step=1/(handles.image_depth-1);
max_step=5*min_step;
set(handles.slider1, ...
    'Enable','on', ...
    'Min',1, ...
    'Max',handles.image_depth, ...
    'Value',1, ...
    'SliderStep', [min_step max_step]);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes behavior_correct_speed wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = behavior_correct_speed_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in all_frames.
% This button is to update all frames following
function all_frames_Callback(hObject, eventdata, handles)
% hObject    handle to all_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pointofchange = handles.framenumber;
handles.velant(pointofchange:end) = ...
    -handles.velant(pointofchange:end);
handles.velpos(pointofchange:end) = ...
    -handles.velpos(pointofchange:end);
handles.centerline(:,2*pointofchange-1:end) = ...
    flipud(handles.centerline(:,2*pointofchange-1:end));

% Update velocity figure
axes(handles.Velocity);
cla;
% plot(handles.velant,'g');
% hold on;
plot(handles.velpos,'k'); hold on;
plot([1 handles.image_depth], [0 0], 'r');
ylim([-handles.mag*handles.med handles.mag*handles.med]);
plot([handles.framenumber handles.framenumber], ylim, 'm');

% Update recording figure
axes(handles.Recording);
cla;
img=imagesc(handles.img_stack{handles.framenumber+handles.istart-1,1});
colormap(gray); hold on;
i=handles.framenumber+handles.istart-1;
plot(handles.dorsal{2*i-1,1}, handles.dorsal{2*i,1}, 'r');
plot(handles.ventral{2*i-1,1}, handles.ventral{2*i,1}, 'r');
plot(handles.centerline(:, 2*i-1), handles.centerline(:, 2*i), 'w', 'linewidth', 1.5);
plot(handles.centerline(1,2*i-1), handles.centerline(1, 2*i), 'og', 'markersize', 4, 'markerfacecolor', 'g');
plot(handles.centerline(end/2,2*i-1), handles.centerline(end/2, 2*i), 'ok', 'markersize', 4, 'markerfacecolor', 'k');
plot(handles.centerline(end,2*i-1), handles.centerline(end, 2*i), 'oy', 'markersize', 4, 'markerfacecolor', 'y');
text(10,10,['Frame #' num2str(handles.framenumber)],'color', 'k');
set(handles.Recording, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0], 'Units', 'pixels', 'DataAspectRatio', [1 1 1]);
set(img,'ButtonDownFcn', 'proof_reading(''ButtonDown_Callback'',gcbo,[],guidata(gcbo))');

beep;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in one_frame.
% This button is to update only the current frame
function one_frame_Callback(hObject, eventdata, handles)
% hObject    handle to one_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pointofchange = handles.framenumber;
handles.velant(pointofchange) = ...
    -handles.velant(pointofchange);
handles.velpos(pointofchange) = ...
    -handles.velpos(pointofchange);
handles.centerline(:,(2*pointofchange-1):2*pointofchange) = ...
    flipud(handles.centerline(:,(2*pointofchange-1):2*pointofchange));


% Update velocity
axes(handles.Velocity);
cla;
% plot(handles.velant,'g');
% hold on;
plot(handles.velpos,'k'); hold on;
plot([1 handles.image_depth], [0 0], 'r');
ylim([-handles.mag*handles.med handles.mag*handles.med]);
plot([handles.framenumber handles.framenumber], ylim, 'm');

% Update recording
axes(handles.Recording);
cla;
img=imagesc(handles.img_stack{handles.framenumber+handles.istart-1,1});
colormap(gray); hold on;
i=handles.framenumber+handles.istart-1;
plot(handles.dorsal{2*i-1,1}, handles.dorsal{2*i,1}, 'r');
plot(handles.ventral{2*i-1,1}, handles.ventral{2*i,1}, 'r');
plot(handles.centerline(:, 2*i-1), handles.centerline(:, 2*i), 'w', 'linewidth', 1.5);
plot(handles.centerline(1,2*i-1), handles.centerline(1, 2*i), 'og', 'markersize', 4, 'markerfacecolor', 'g');
plot(handles.centerline(end/2,2*i-1), handles.centerline(end/2, 2*i), 'ok', 'markersize', 4, 'markerfacecolor', 'k');
plot(handles.centerline(end,2*i-1), handles.centerline(end, 2*i), 'oy', 'markersize', 4, 'markerfacecolor', 'y');
text(10,10,['Frame #' num2str(handles.framenumber)],'color', 'k');
set(handles.Recording, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0], 'Units', 'pixels', 'DataAspectRatio', [1 1 1]);
set(img,'ButtonDownFcn', 'proof_reading(''ButtonDown_Callback'',gcbo,[],guidata(gcbo))');

beep;

% Update handles structure
guidata(hObject, handles);

% function behavior_correct_speed_WindowKeyPressFcn(hObject, eventdata, handles)
%  % determine the key that was pressed 
%  keyPressed = eventdata.Key;
%  if strcmpi(keyPressed,'x')
%      % set focus to the button
%      uicontrol(handles.one_frame);
%      % call the callback
%      pushbutton1_Callback(handles.one_frame,[],handles);
%  end
 
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.framenumber=round(get(hObject,'Value'));
% set(handles.text1, 'String',strcat(num2str(handles.frame_number+handles.istart-1),'/',num2str(handles.sections),'(',handles.filename,')'));

% Update velocity
axes(handles.Velocity); 
cla;
% plot(handles.velant,'g');
% hold on;
plot(handles.velpos,'k'); hold on;
plot([1 handles.image_depth], [0 0], 'r');
ylim([-handles.mag*handles.med handles.mag*handles.med]);
plot([handles.framenumber handles.framenumber], ylim, 'm');

% Update recording
axes(handles.Recording);
cla;
img=imagesc(handles.img_stack{handles.framenumber+handles.istart-1,1});
colormap(gray); hold on;
i=handles.framenumber+handles.istart-1;
plot(handles.dorsal{2*i-1,1}, handles.dorsal{2*i,1}, 'r');
plot(handles.ventral{2*i-1,1}, handles.ventral{2*i,1}, 'r');
plot(handles.centerline(:, 2*i-1), handles.centerline(:, 2*i), 'w', 'linewidth', 1.5);
plot(handles.centerline(1,2*i-1), handles.centerline(1, 2*i), 'og', 'markersize', 4, 'markerfacecolor', 'g');
plot(handles.centerline(end/2,2*i-1), handles.centerline(end/2, 2*i), 'ok', 'markersize', 4, 'markerfacecolor', 'k');
plot(handles.centerline(end,2*i-1), handles.centerline(end, 2*i), 'oy', 'markersize', 4, 'markerfacecolor', 'y');
text(10,10,['Frame #' num2str(handles.framenumber)],'color', 'k');
set(handles.Recording, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0], 'Units', 'pixels', 'DataAspectRatio', [1 1 1]);
set(img,'ButtonDownFcn', 'proof_reading(''ButtonDown_Callback'',gcbo,[],guidata(gcbo))');

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


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


% --- Executes on button press in export_data.
function export_data_Callback(hObject, eventdata, handles)
% hObject    handle to export_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

assignin('base','velocity_corrected_anterior',handles.velant);
assignin('base','velocity_corrected_posterior',handles.velpos);
assignin('base','centerline_corrected',handles.centerline);

fprintf(['corrected speed is exported for ' handles.filename '.\n']);
fprintf('speed manual curation is done. \n');
