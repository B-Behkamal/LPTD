function varargout = LPTD_GUI_2(varargin)


% LPTD_GUI_2 MATLAB code for LPTD_GUI_2.fig
%      LPTD_GUI_2, by itself, creates a new LPTD_GUI_2 or raises the existing
%      singleton*.
%
%      H = LPTD_GUI_2 returns the handle to a new LPTD_GUI_2 or the handle to
%      the existing singleton*.
%
%      LPTD_GUI_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LPTD_GUI_2.M with the given input arguments.
%
%      LPTD_GUI_2('Property','Value',...) creates a new LPTD_GUI_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LPTD_GUI_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LPTD_GUI_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LPTD_GUI_2

% Last Modified by GUIDE v2.5 30-Jan-2022 09:33:17

% Begin initialization code - DO NOT EDIT


gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LPTD_GUI_2_OpeningFcn, ...
    'gui_OutputFcn',  @LPTD_GUI_2_OutputFcn, ...
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


% --- Executes just before LPTD_GUI_2 is made visible.

function LPTD_GUI_2_OpeningFcn(hObject, eventdata, handles, varargin)

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

% varargin   command line arguments to LPTD_GUI_2 (see VARARGIN)

% Choose default command line output for LPTD_GUI_2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LPTD_GUI_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);



radiobutton_helix = findobj(0, 'tag', 'radiobutton_helix');
radiobutton_sheet = findobj(0, 'tag', 'radiobutton_sheet');

uipanel1 = findobj(0, 'tag', 'uipanel1');
uipanel2 = findobj(0, 'tag', 'uipanel2');

if get(radiobutton_helix,'Value')==1
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'off')
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'on')
    
    
elseif get(radiobutton_sheet,'Value')==1
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'off')
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'on')
    
end

uipanel1 = findobj(0, 'tag', 'uipanel1');
set(findall(uipanel1, '-property', 'enable'), 'enable', 'off')

uipanel2 = findobj(0, 'tag', 'uipanel2');
set(findall(uipanel2, '-property', 'enable'), 'enable', 'off')



% --- Outputs from this function are returned to the command line.
function varargout = LPTD_GUI_2_OutputFcn(hObject, eventdata, handles)

% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;



function et_pdb_strand_Callback(hObject, eventdata, handles)

% hObject    handle to et_pdb_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_pdb_strand as text
%        str2double(get(hObject,'String')) returns contents of et_pdb_strand as a double


% --- Executes during object creation, after setting all properties.

function et_pdb_strand_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_pdb_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_stick_strand_Callback(hObject, eventdata, handles)

% hObject    handle to et_stick_strand (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_stick_strand as text
%        str2double(get(hObject,'String')) returns contents of et_stick_strand as a double


% --- Executes during object creation, after setting all properties.

function et_stick_strand_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_stick_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_chain_strand_Callback(hObject, eventdata, handles)

% hObject    handle to et_chain_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_chain_strand as text
%        str2double(get(hObject,'String')) returns contents of et_chain_strand as a double


% --- Executes during object creation, after setting all properties.

function et_chain_strand_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_chain_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.

function pushbutton2_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

et_pdb_strand = findobj(0, 'tag', 'et_pdb_strand');
pdb_id=get(et_pdb_strand,'String');
pdb_id=string(pdb_id);
protein_data.name=pdb_id;

et_stick_strand = findobj(0, 'tag', 'et_stick_strand');
stick_strands=get(et_stick_strand,'String');
stick_strands=string(stick_strands);
protein_data.stick_strands=stick_strands;

et_chain_strand = findobj(0, 'tag', 'et_chain_strand');
chain=get(et_chain_strand,'String');
chain='A'
chain=string(chain);
protein_data.chain=chain;


radiobutton_helix = findobj(0, 'tag', 'radiobutton_helix');
helix_info=get(radiobutton_helix,'value');
helix_info=double(helix_info(1));
protein_data.helix_info=helix_info;


radiobutton_sheet = findobj(0, 'tag', 'radiobutton_sheet');
sheet_info=get(radiobutton_sheet,'value');
sheet_info=double(sheet_info(1));
protein_data.sheet_info=sheet_info;

if sheet_info==1
    final_topology=LPTD_Function_2(protein_data,helix_info,sheet_info)
    final_topology_sheet= final_topology;
    x=[[final_topology_sheet(:).num_strand]',[final_topology_sheet(:).num_stick]',[final_topology_sheet(:).Direction]'];
    m1=num2cell(x);
    for i=1:size(m1,1)
        if m1{i,3}==0
            m1{i,3}=' ';
        end
    end
    
    
    Topology_strand = findobj(0, 'tag', 'Topology_strand');
    set(Topology_strand,'Data',m1)
    
end


function et_pdb_helix_Callback(hObject, eventdata, handles)

% hObject    handle to et_pdb_helix (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_pdb_helix as text
%        str2double(get(hObject,'String')) returns contents of et_pdb_helix as a double


% --- Executes during object creation, after setting all properties.
function et_pdb_helix_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_pdb_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_stick_helix_Callback(hObject, eventdata, handles)

% hObject    handle to et_stick_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_stick_helix as text
%        str2double(get(hObject,'String')) returns contents of et_stick_helix as a double


% --- Executes during object creation, after setting all properties.
function et_stick_helix_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_stick_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_chain_helix_Callback(hObject, eventdata, handles)

% hObject    handle to et_chain_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_chain_helix as text
%        str2double(get(hObject,'String')) returns contents of et_chain_helix as a double


% --- Executes during object creation, after setting all properties.
function et_chain_helix_CreateFcn(hObject, eventdata, handles)

% hObject    handle to et_chain_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

check_helix = findobj(0, 'tag', 'radiobutton_helix');
check_sheet = findobj(0, 'tag', 'radiobutton_sheet');



if get(check_helix, 'Value') == 1
    et_pdb = findobj(0, 'tag', 'et_pdb_helix');
    pdb_id=get(et_pdb,'String');
    pdb_id=string(pdb_id);
    
    protein_data.name=pdb_id;
    
    et_stick_helix = findobj(0, 'tag', 'et_stick_helix');
    stick_hlces=get(et_stick_helix,'String');
    stick_hlces=string(stick_hlces);
    protein_data.stick_hlces=stick_hlces;
    
    et_chain = findobj(0, 'tag', 'et_chain_helix');
    chain=get(et_chain,'String');
    chain='A';
    chain=string(chain);
    protein_data.chain=chain;
    
    
    radiobutton_helix = findobj(0, 'tag', 'radiobutton_helix');
    helix_info=get(radiobutton_helix,'value');
    helix_info=double(helix_info(1));
    protein_data.helix_info=helix_info;
    
    sheet_tag = findobj(0,'tag', 'radiobutton_sheet');
    sheet_info = double(get(sheet_tag, 'Value'));
end


if get(check_sheet, 'Value') == 1
    et_pdb = findobj(0, 'tag', 'et_pdb_strand');
    pdb_id=get(et_pdb,'String');
    pdb_id=string(pdb_id);
    protein_data.name=pdb_id;
    
    et_chain = findobj(0, 'tag', 'et_chain_strand');
    chain=get(et_chain,'String');
    chain='A';
    chain=string(chain);
    protein_data.chain=chain;
    
    
    et_stick_strand = findobj(0, 'tag', 'et_stick_strand');
    stick_strands=get(et_stick_strand,'String');
    stick_strands=string(stick_strands);
    protein_data.stick_strands=stick_strands;
    
    
    radiobutton_sheet = findobj(0, 'tag', 'radiobutton_sheet');
    sheet_info=get(radiobutton_sheet,'value');
    sheet_info=double(sheet_info(1));
    protein_data.sheet_info=sheet_info;
    
    helix_tag = findobj(0,'tag', 'radiobutton_helix');
    helix_info = double(get(helix_tag, 'Value'));
end


[final_topology]=LPTD_Function_2(protein_data,helix_info,sheet_info);

if helix_info == 1 && sheet_info == 0
    final_topology_helix = final_topology;
    
    x_helix=[[final_topology_helix(:).num_helix]',[final_topology_helix(:).num_stick]',[final_topology_helix(:).Direction]'];
    m_helix=num2cell(x_helix);
    for i=1:length(m_helix)
        if m_helix{i,3}==0
            m_helix{i,3}=' ';
        end
    end

    
    Topology_helix = findobj(0, 'tag', 'Topology_helix');
    set(Topology_helix,'Data',m_helix)
    
elseif sheet_info ==1 && helix_info == 0
    
    final_topology_sheet = final_topology;
    
    x=[[final_topology_sheet(:).num_strand]',[final_topology_sheet(:).num_stick]',[final_topology_sheet(:).Direction]'];
    m_sheet=num2cell(x);
    
    for i=1:size(m_sheet,1)
        if m_sheet{i,3}==0
            m_sheet{i,3}=' ';
        end
    end
    
    Topology_strand = findobj(0, 'tag', 'Topology_strand');
    set(Topology_strand,'Data',m_sheet)
    
elseif helix_info == 1 && sheet_info == 1
    
    
    final_topology_helix = cell2mat(final_topology(1));
    x_helix=[[final_topology_helix(:).num_helix]',[final_topology_helix(:).num_stick]',[final_topology_helix(:).Direction]'];
    m_helix=num2cell(x_helix);
    
    for i=1:length(m_helix)
        if m_helix{i,3}==0
            m_helix{i,3}=' ';
        end
    end
    
    Topology_helix = findobj(0, 'tag', 'Topology_helix');
    set(Topology_helix,'Data',m_helix)
    
    
    final_topology_sheet = cell2mat(final_topology(2));
    x=[[final_topology_sheet(:).num_strand]',[final_topology_sheet(:).num_stick]',[final_topology_sheet(:).Direction]'];
    m_sheet=num2cell(x);
    
    for i=1:size(m_sheet,1)
        if m_sheet{i,3}==0
            m_sheet{i,3}=' ';
        end
    end
    
    Topology_strand = findobj(0, 'tag', 'Topology_strand');
    set(Topology_strand,'Data',m_sheet)
    
    
end


Topology_helix = findobj(0, 'tag', 'Topology_helix');
set(Topology_helix,'ColumnName',{'#Helix', '#Stick', '#Direction'})

Topology_strand = findobj(0, 'tag', 'Topology_strand');
set(Topology_strand,'ColumnName',{'#Strand', '#Stick', '#Direction'})



% --- Executes during object creation, after setting all properties.
function chechboxes_CreateFcn(hObject, eventdata, handles)

% hObject    handle to chechboxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in chechboxes.
function chechboxes_SelectionChangedFcn(hObject, eventdata, handles)

% hObject    handle to the selected object in chechboxes
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

radiobutton_helix = findobj(0, 'tag', 'radiobutton_helix');
radiobutton_sheet = findobj(0, 'tag', 'radiobutton_sheet');

uipanel1 = findobj(0, 'tag', 'uipanel1');
uipanel2 = findobj(0, 'tag', 'uipanel2');
if get(radiobutton_helix,'Value')==1
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'off')
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'on')
    
    
elseif get(radiobutton_sheet,'Value')==1
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'off')
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'on')
    
end



function et_topology_helix_Callback(hObject, eventdata, handles)

% hObject    handle to et_topology_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_topology_helix as text
%        str2double(get(hObject,'String')) returns contents of et_topology_helix as a double


% --- Executes during object creation, after setting all properties.
function et_topology_helix_CreateFcn(hObject, eventdata, handles)


% hObject    handle to et_topology_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function et_topology_strand_Callback(hObject, eventdata, handles)


% hObject    handle to et_topology_strand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of et_topology_strand as text
%        str2double(get(hObject,'String')) returns contents of et_topology_strand as a double


% --- Executes during object creation, after setting all properties.

function et_topology_strand_CreateFcn(hObject, eventdata, handles)


% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in Topology_helix.
function Topology_helix_CellEditCallback(hObject, eventdata, handles)


% hObject    handle to Topology_helix (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Topology_helix_CreateFcn(hObject, eventdata, handles)


% hObject    handle to Topology_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in showresults.
function showresults_Callback(hObject, eventdata, handles)


% hObject    handle to showresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobutton_helix.
function radiobutton_helix_Callback(hObject, eventdata, handles)


% hObject    handle to radiobutton_helix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_helix

checkbox_helix = findobj(0, 'tag', 'radiobutton_helix');


check_helix = get(checkbox_helix, 'Value');

uipanel1 = findobj(0, 'tag', 'uipanel1');

if check_helix == 1
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'on')
else
    set(findall(uipanel1, '-property', 'enable'), 'enable', 'off')
end


% --- Executes on button press in radiobutton_sheet.
function radiobutton_sheet_Callback(hObject, eventdata, handles)


% hObject    handle to radiobutton_sheet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_sheet


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over radiobutton_helix.

checkbox_sheet = findobj(0, 'tag', 'radiobutton_sheet');

check_sheet = get(checkbox_sheet, 'Value');

uipanel2 = findobj(0, 'tag', 'uipanel2');

if check_sheet == 1
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'on')
else
    set(findall(uipanel2, '-property', 'enable'), 'enable', 'off')
end


% --------------------------------------------------------------------
function chechboxes_ButtonDownFcn(hObject, eventdata, handles)

% hObject    handle to chechboxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function precision_Callback(hObject, eventdata, handles)

% hObject    handle to precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of precision as text
%        str2double(get(hObject,'String')) returns contents of precision as a double


% --- Executes during object creation, after setting all properties.
function precision_CreateFcn(hObject, eventdata, handles)

% hObject    handle to precision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function recall_Callback(hObject, eventdata, handles)

% hObject    handle to recall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recall as text
%        str2double(get(hObject,'String')) returns contents of recall as a double


% --- Executes during object creation, after setting all properties.

function recall_CreateFcn(hObject, eventdata, handles)

% hObject    handle to recall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fmeasure_Callback(hObject, eventdata, handles)

% hObject    handle to fmeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fmeasure as text
%        str2double(get(hObject,'String')) returns contents of fmeasure as a double



% --- Executes during object creation, after setting all properties.

function fmeasure_CreateFcn(hObject, eventdata, handles)

% hObject    handle to fmeasure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function accuracy_Callback(hObject, eventdata, handles)
% hObject    handle to accuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of accuracy as text
%        str2double(get(hObject,'String')) returns contents of accuracy as a double



% --- Executes during object creation, after setting all properties.

function accuracy_CreateFcn(hObject, eventdata, handles)

% hObject    handle to accuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rank_Callback(hObject, eventdata, handles)
% hObject    handle to rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rank as text
%        str2double(get(hObject,'String')) returns contents of rank as a double



% --- Executes during object creation, after setting all properties.

function rank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
