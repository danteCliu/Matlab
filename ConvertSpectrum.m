function varargout = ConvertSpectrum(varargin)
% CONVERTSPECTRUM MATLAB code for ConvertSpectrum.fig
%basic assumes :
 % 0, sample and condition is exactly the same, difference is incident wave
 % 1, Lab XRD Intensity is counts per seconds
 % 2, APS spectrum intensity is count per seconds
 % 3, use powder diffraction intensity formula : I = I0 (lambda^3*e^4*V/32pi
 % Rm^2c^4vo^2)*pF^2LpA(theta)e^(-2M). exclude factor that from the
 % sample and the experiment setup.
 % 4, bragg law used 2dsin(theta) = n*lambda
 % 5, from KeV to Lambda (AEDIT) use relation: 1KeV = 1.23984 (nm) = 12.398425 AEDIT
 % E = hc/¦Ë so wavelenth is reciprocal to lambda(AEDIT) = 12.398425/keV
 % ¦Í = aEdit/¦Ë
 % h = Planck's constant = 4.135667516 x 10-15 eV*s
 % aEdit = speed of light = 299792458 m/s
 %
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConvertSpectrum

% Last Modified by GUIDE v2.5 25-Mar-2016 09:59:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConvertSpectrum_OpeningFcn, ...
                   'gui_OutputFcn',  @ConvertSpectrum_OutputFcn, ...
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


% --- Executes just before ConvertSpectrum is made visible.
function ConvertSpectrum_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConvertSpectrum (see VARARGIN)
global udata;
% Choose default command line output for ConvertSpectrum
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
udata.wavefactor = 12.398425; %KeV to Lambda (Angstrom)
udata.Lablambda = 1.544; %
% LogoAxes = handles.LogoAxes;
% I = imread('logo.jpg');
% imshow(I,'parent',LogoAxes);

% UIWAIT makes ConvertSpectrum wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ConvertSpectrum_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in LoadXRDfile.
function LoadXRDfile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadXRDfile (see GCBO)
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata;
% read data
[filename, pathname] = uigetfile('*.*','pick a XRD file');
if ispc
    foldername = regexp(pathname,'\','split');
    foldername = [char(foldername(end-1)),'\',filename];
else
    foldername = regexp(pathname,'/','split');
    foldername = [char(foldername(end-1)),'/',filename];
end

set(handles.XRDfile,'string',foldername);
cd (pathname);

fid = fopen(filename);

 s1 = fscanf(fid,'%s',(1));
 s2 = strfind(s1,';'); %find ';' in scaned string 
 s3 = strfind(s1,','); % if a csv file
 format long
if isempty(s2)&&isempty(s3)
    
    fseek(fid,0,-1); %get to the start point
    % skip lines to be added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xrd = fscanf(fid,'%f %f',[2 inf]);
    
 elseif ~isempty(s3)
    
    fseek(fid,0,-1);
    
    xrd = fscanf(fid,'%f,%f',[2 inf]);
 else
    fseek(fid,0,-1);
    
    xrd = fscanf(fid,'%f;%f',[2 inf]);
    
end

fclose(fid);
udata.xrd = xrd;

% --- Executes on button press in Load_spectrum.
function Load_spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to Load_spectrum (see GCBO)
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global udata;
% read data
[filename, pathname] = uigetfile('*.txt','pick a spectrum file');
if ispc
    foldername = regexp(pathname,'\','split');
    foldername = [char(foldername(end-1)),'\',filename];
else
    foldername = regexp(pathname,'/','split');
    foldername = [char(foldername(end-1)),'/',filename];
end

 set(handles.specfile,'string',foldername);

 cd (pathname);

fid = fopen(filename);
 spectrum = fscanf(fid,'%f %f',[2 inf]);
fclose(fid);
udata.spec = spectrum;
% --- Executes on selection change in popupmenu1.


% --- Executes on button press in Polt_Lab.
function Polt_Lab_Callback(hObject, eventdata, handles)
% hObject    handle to Polt_Lab (see GCBO)
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata;
xrd = udata.xrd;  % 2*inf array to plot

 xaxe = xrd( 1 , :);
 yaxe = xrd( 2 , :);
 figure('Tag','XRDprofile');
 plot(xaxe, yaxe,':');
   title('Lab XRD profile','FontSize',20); 
   xlabel('2\theta','FontSize',20);
   ylabel('Intensity','FontSize',20);
   
function minheightEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minheightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata
% Hints: get(hObject,'String') returns contents of minheightEdit as text
%        str2double(get(hObject,'String')) returns contents of minheightEdit as a double
minh = get(handles.minheightEdit,'string');
udata.minh = str2double(minh);

% --- Executes during object creation, after setting all properties.
function minheightEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minheightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mindistanceEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mindistanceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata
% Hints: get(hObject,'String') returns contents of mindistanceEdit as text
%        str2double(get(hObject,'String')) returns contents of mindistanceEdit as a double
mind = get(handles.mindistanceEdit,'string');
udata.mind = str2double(mind);

% --- Executes during object creation, after setting all properties.
function mindistanceEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mindistanceEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in Polt_Spectrum.
function Polt_Spectrum_Callback(hObject, eventdata, handles)
% hObject    handle to Polt_Spectrum (see GCBO)
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata;
spec = udata.spec;  % 2*inf array to plot

xaxe = spec( 1 , :);
yaxe = spec( 2 , :);
 figure('Tag','xrd_spectrum');
 plot(xaxe, yaxe);
   title('Spectrum','FontSize',20); 
   xlabel('Energy (keV)','FontSize',20);
   ylabel('Flux','FontSize',20);



 
% --- Executes on button press in Convert.
function Convert_Callback(hObject, eventdata, handles)
% hObject    handle to Convert (see GCBO)
% eventdata  reserved - to be defined in aEdit future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global udata;
  xrd  = udata.xrd;     %data from lab [2 inf]
  spec = udata.spec;    %wave difference [2 inf]  in column order
  wafac = udata.wavefactor;
  lambda0 = udata.Lablambda;
%   Inci_inten = udata.LabInci;
  minh = udata.minh;
  mind = udata.mind;
  % cacaulate sample const from lab data 
  theta_a = xrd(1, : );
  theta_b = theta_a./2 ; % bragg angle theta in data.
  inten = xrd(2 , : );  % intensity in data.
  
 %    Lab_K = inten./(Inci_inten*lambda0.^3);  %Intensity formula k about samples 
   [~, inde] = findpeaks(inten,'MINPEAKHEIGHT',minh,'MINPEAKDISTANCE',mind);
     % plot peaks found out
     peakx = theta_a(inde); 
     peaky = inten(inde);
     figure('Tag','Xrdpeaks');
     plot(theta_a , inten,':',peakx,peaky,'ro');
     title('Peaks found out');
     xlabel('2\theta');
     ylabel('Intensity');
    % pick theta strengthen' d-spacing theta
     theta_d = theta_b(inde);
   
     % cacaulate d spacing assume lab monochrome xray in the peaks
     d_lab = lambda0./(2*sind(theta_d)); 
        
  % convert KeV to Angstrom
    flux = spec(2, :);
    energy_a = spec(1, : );   % to be part
    spec_lamb = wafac./energy_a;
    figure('Tag','lamdatoflux');
    plot(spec_lamb,flux);
    xlabel('lambda(ang)');
    ylabel('flux');
    %read spectrum wavelenth and get theta value of each d-spaceing
    [rown m] = size(d_lab);         %1*m array.
    if rown == 1
        disp('your spectrum is right!');
    else
        disp('wrong in read of spectrum!');
    end
    % set all sin theta value of every d spaceing
    [~,speccol] = size(spec_lamb);
    theta_sindata = zeros(m,speccol);
    for ii = 1:1:m
          
        
        theta_sin = spec_lamb./(2*d_lab(ii));       % theta from small to big
        judge = theta_sin <= ones(size(theta_sin)); % bragg angle theta 
        theta_sinR = theta_sin.*judge;
        theta_sindata(ii,:) = theta_sinR;
    end
    % collect all point
    %     [nn mm] = size(xrd); %2*mm arrays
    %     pres = (max(theta_a)-min(theta_a))/mm;
        dataout = cell(1,m);
    for ii = 1:1:m
          
        
        con_theta = 2*asind(theta_sindata(ii,:)); % theta to 2theta
        con_thetarand = roundn(con_theta,-2);     %presicion
        Inten = flux.*spec_lamb.^3;               % 16/3/15consider theta involve
        
        profile1 = [con_thetarand;Inten]; % the profile array
         
        dataout{1,ii} = profile1; 
    end
    %% find out same theta accumulate intensity
     % creat a format theta and compare theta
     endang = 90;
     rown = ((90-0)/0.01) + 1;
     assume_theta = (0:0.01:endang);
     data_all_inten = zeros(m,rown);
     out_d = zeros(1,9001);
    for ii = 1:1:m
        data_theta = dataout{1,ii}(1,:);
        data_inten = dataout{1,ii}(2,:);
        same_theta = intersect(data_theta,assume_theta);
        [~,mm] = size(same_theta);
        
        for jj = 1:1:mm
            
            index_the = find(assume_theta == same_theta(jj));% as data_theta
            index_inten = find(data_theta == same_theta(jj));
            
            
            out_d(index_the) = data_inten(index_inten(1));
        end
        
        data_all_inten(ii,:) = out_d;
    end
    converted_X = assume_theta';
    converted_Y = sum(data_all_inten)';
    XY = [converted_X converted_Y];    % 2*inf array
    figure('Tag','final plot');
    plot(converted_X,converted_Y,':');
    title('simulated profile');
    xlabel('2theta');
    ylabel('Intensity');
    save('output_file.dat','-ascii','XY');
    %% end of code not include temperature factor absorption factor. mostally 2dsin(theta) = n lamb n not in consider
        % and d-spaceing just count from the lab xrd profile peaks, do not counts all possibilties should miss some peaks
        % in some case, to be fixed with .cif file cout as many d-spaceing as
        % possible.
        
  


% --- Executes on button press in closefig.
function closefig_Callback(hObject, eventdata, handles)
% hObject    handle to closefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
delete(allPlots);
