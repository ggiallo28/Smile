function varargout = ROIDraw(varargin)
% ROIDraw M-file for ROIDraw.fig
%      ROIDraw, by itself, creates a new ROIDraw or raises the existing
%      singleton*.
%
%      H = ROIDraw returns the handle to a new ROIDraw or the handle to
%      the existing singleton*.
%
%      ROIDraw('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROIDraw.M with the given input arguments.
%
%      ROIDraw('Property','Value',...) creates a new ROIDraw or
%      raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROIDraw_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROIDraw_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROIDraw

% Last Modified by GUIDE v2.5 25-Nov-2008 11:45:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROIDraw_OpeningFcn, ...
                   'gui_OutputFcn',  @ROIDraw_OutputFcn, ...
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


% --- Executes just before ROIDraw is made visible.
function ROIDraw_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROIDraw (see VARARGIN)

% Choose default command line output for ROIDraw
handles.output = hObject;
%define variables
axis off;
handles.a_x=0;                      % coordinates of four corners of rectangle
handles.a_y=0;                  
handles.b_x=0;
handles.b_y=0;
handles.c_x=0;
handles.c_y=0;
handles.d_x=0;
handles.d_y=0;
handles.X=0;                        % array which stores all X's of corners
handles.Y=0;                        % array which stores all Y's of corners
handles.filename = 'blank';         % stores DICOM filename
handles.output = hObject;
handles.atob_theta_indegrees = 0;   % angle of a to b, measured in degrees
handles.btoc_theta_indegrees = 0;   % angle of b to c, measured in degrees
handles.is_rectangle = 0;           % boolean, true if rectangle is vertical
handles.x_min = 0;                  
handles.y_min = 0;
handles.x_max = 0;
handles.y_max = 0;        
handles.roi_defined =0;             % boolean, true if roi has been drawn
handles.info = '';                  % to store dicom info, not implemented yet
handles.preview = 0;                % stores handle to preview figure
handles.imagecopy='';               
handles.ones='';                    % stores array of ones to be used to crop image
handles.just_filename='';           % file name without path name
handles.mat_file='';                % .mat file loaded for colocalized ROI
handles.original = '';              % original image
handles.use_load = 0;               
handles.roi = '';
handles.mask2= '';
handles.line = '';                  % line used to draw rectangle
handles.file_open = 0;              % boolean, true if file has been opened
handles.x_vals = '';
handles.y_vals = '';
handles.was_previewed = 0;
handles.save_directory = cd;
handles.workspace = cd;
handles.opendir = cd;
handles.loaddir = cd;
handles.zoom='';
handles.current_patient='';
handles.userset=0;
handles.figure_image='';
handles.height=0;
handles.width=0;

save_text = strcat(handles.save_directory,'\');
set(handles.save_directory_text,'String',save_text);


% results
handles.result = '';                % resulting cropped ROI image
handles.mask = '';                  % resulting mask image
handles.image = '';                 % original image

set(handles.status_text,'String','Open a file to begin selecting an ROI.');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ROIDraw wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ROIDraw_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in draw_roi.
function draw_roi_Callback(hObject, eventdata, handles)
% hObject    handle to draw_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.is_rectangle = 0;
if(~handles.file_open)
    errordlg('Please open a file before attempting to draw an ROI.');
else % file has been opened
    if(handles.roi_defined)
        % reset everything if an roi was previously defined
        v = caxis;
        r = axis;

        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);
    end
    set(handles.status_text,'String','Draw a rectangle with four points. Hit Enter once you have finished.');
    
    % get the lines drawn by the user
    line_x=[];
    line_y=[];
    bad_rectangle = 1;
    while(bad_rectangle == 1)
        [line_x,line_y] = getline();
        if (size(line_x) < 3)
            e = errordlg('You must define at least 3 points to create your rectangle. Try again.');
            waitfor(e);
        else
            bad_rectangle=0;
        end        
    end
    
    handles.a_x = line_x(1);
    handles.a_y = line_y(1);
    handles.b_x = line_x(2);
    handles.b_y = line_y(2);
    tempc_x = line_x(3);
    tempc_y = line_y(3);
    
    % find the difference between point b and point a 
    % so that we don't draw a parallelogram
    diff2_x = handles.b_x - handles.a_x;
    diff2_y = handles.b_y - handles.a_y;
     
    sign2_x=0;
    sign2_y=0;
    % we only care about its sign
    if(diff2_x ~= 0)
        sign2_x = diff2_x/(sqrt(diff2_x^2));
    end
    if(diff2_y ~= 0)
        sign2_y = diff2_y/(sqrt(diff2_y^2));
    end

    % find the slopes of our lines
    m_atob = (handles.b_y - handles.a_y)/(handles.b_x - handles.a_x);
    m_btoc = -1/m_atob;
    b_atob = handles.b_y - m_atob*handles.b_x;
    
    % determine which side of segment atob point c lies
    lies_above=0;
    if(tempc_y >= (m_atob*tempc_x + b_atob))
        lies_above=1;
    end
      
    % find the angle of our lines
    atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
    handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;

    btoc_theta_inradians = atan(abs(m_btoc));
    handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;

    % find the lengths of our rectangle
    %length_atob = sqrt((handles.a_y-handles.b_y)^2 + (handles.a_x - handles.b_x)^2);
    length_btoc = sqrt((handles.b_y-tempc_y)^2 + (handles.b_x-tempc_x)^2);

    % we already have points and b, but now we must calculate c and d
    % based on the slope of segment ab and on where point c lies in
    % relation to segment ab    
    if(m_atob >0)
        if(lies_above)            
            handles.c_x = handles.b_x - abs(length_btoc * cos(btoc_theta_inradians));
            handles.c_y = handles.b_y + abs(length_btoc * sin(btoc_theta_inradians));
            handles.d_x = handles.a_x - abs(length_btoc * cos(btoc_theta_inradians));
            handles.d_y = handles.a_y + abs(length_btoc * sin(btoc_theta_inradians));
        else
            handles.c_x = handles.b_x + abs(length_btoc * cos(btoc_theta_inradians));
            handles.c_y = handles.b_y - abs(length_btoc * sin(btoc_theta_inradians));
            handles.d_x = handles.a_x + abs(length_btoc * cos(btoc_theta_inradians));
            handles.d_y = handles.a_y - abs(length_btoc * sin(btoc_theta_inradians));
        end           
    else
        if(lies_above)            
            handles.c_x = handles.b_x + abs(length_btoc * cos(btoc_theta_inradians));
            handles.c_y = handles.b_y + abs(length_btoc * sin(btoc_theta_inradians));
            handles.d_x = handles.a_x + abs(length_btoc * cos(btoc_theta_inradians));
            handles.d_y = handles.a_y + abs(length_btoc * sin(btoc_theta_inradians));
        else
            handles.c_x = handles.b_x - abs(length_btoc * cos(btoc_theta_inradians));
            handles.c_y = handles.b_y - abs(length_btoc * sin(btoc_theta_inradians));
            handles.d_x = handles.a_x - abs(length_btoc * cos(btoc_theta_inradians));
            handles.d_y = handles.a_y - abs(length_btoc * sin(btoc_theta_inradians));
        end         
    end

    % special case, handle perpendicular rectangles
    if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
        % we have a normal non-rotated rectangle
        handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
        handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
        handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
        handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

        % draw the rectangle
        rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
        handles.is_rectangle = 1;

        handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
        handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
        handles.mask = roipoly(handles.image, handles.X, handles.Y);
    else % the rectangle is not a vertical rectangle, has an angle
        handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
        handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

        % draw the rectangle using the points we've found
        line(handles.X,handles.Y,'Marker','.','LineStyle','-');
        handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
        handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];
        handles.mask = roipoly(handles.image, handles.x_vals, handles.y_vals);
        
    end

    handles.roi_defined =1;
    set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file before attempting to reset.');
else % the file is open, reset everything and re-display
    
    v = caxis;
    r = axis;
        
    handles.image = dicomread(handles.filename);
    handles.dim=size(handles.image);
    handles.cropped=handles.image; % let's make a copy
    handles.result=handles.image;
    handles.figure_image=imshow(handles.image,'DisplayRange',[]);
    handles.output = hObject;
    handles.was_previewed=0;
    handles.roi_defined=0;
    set(handles.status_text,'String','Click Draw ROI to select an ROI region, or click Load ROI to colocalize an ROI from another image.');

    caxis([v(1) v(2)]);
    axis([r(1) r(2) r(3) r(4)]);
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
% hObject    handle to preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw or load an ROI before attempting to preview.');
else
    if(handles.roi_defined)
        if(handles.is_rectangle)
            % save height and width
            handles.height = handles.y_max - handles.y_min;
            handles.width = handles.x_max - handles.x_min;

            % adjust sampling box determined by pixelspacing
            spacing = handles.info.PixelSpacing;
            spacing_factorx = spacing(2)/.5;
            spacing_factory = spacing(1)/.5;
            pseudo_height = handles.height*spacing_factory;
            pseudo_width = handles.width*spacing_factorx;

            % determine the number of pixels in our resulting ROI 
            % by rounding the width and length
            result = ones(round(pseudo_height),round(pseudo_width));
            translate = [handles.y_min;handles.x_min];

            %initialize rotation data array as emptys
            X=(1:size(handles.image,2));
            Y=(1:size(handles.image,1));
            Z=handles.image;         

            new_pointsx = result;
            new_pointsy = result;
            for i=1:size(result,1),
                for j=1:size(result,2),
                    new_point = [1/spacing_factory*i;1/spacing_factorx*j];
                    new_point = new_point+translate;
                    new_pointsx(i,j) = new_point(2);
                    new_pointsy(i,j) = new_point(1);
                end
            end

            result = interp2(X,Y,Z,new_pointsx,new_pointsy,'bicubic');

            handles.result = result;   
            handles.preview = figure;
            imshow(handles.result,'DisplayRange',[],'InitialMagnification',300);
        else
            if(~handles.was_previewed)                
                % variables used for rotation
                corners_x = handles.x_vals;
                corners_y = handles.y_vals;

                %find the bottom-left corner of the rectangle
                sorted_corners_y = sort(corners_y);
                sorted_corners_x = sort(corners_x);

                % find the lowest y coordinate
                lowest_corner = sorted_corners_y(1);        
                leftmost_corner = sorted_corners_x(4);
                bottom_left_corner = [lowest_corner;leftmost_corner];

                % find the lowest x coordinate that matches with the lowest y coordinate
                for i=1:4;                                  
                    if(corners_y(i) == lowest_corner && corners_x(i) <= leftmost_corner)
                        leftmost_corner = corners_x(i);
                        bottom_left_corner = [lowest_corner;leftmost_corner];
                    end
                end

                % find right-most lowest bottom corner
                rightmost_corner = sorted_corners_x(4);
                lowest_corner = sorted_corners_y(4);
                bottom_right_corner = [lowest_corner;rightmost_corner];
                for i=1:4;
                    if(corners_x(i) == rightmost_corner && corners_y(i) <= lowest_corner)
                        lowest_corner = corners_y(i);
                        bottom_right_corner = [lowest_corner;rightmost_corner];
                    end
                end

                % find left most corner not equal to bottom left corner
                leftmost_corner = sorted_corners_x(1);
                upper_left_corner = [0;leftmost_corner];
                for i=1:4;
                    if(corners_x(i) == leftmost_corner && corners_y(i) ~= bottom_left_corner(1))
                        upper_left_corner = [corners_y(i);leftmost_corner];
                    end
                end

                % set rotation angle
                toright_angle = atand(abs((bottom_left_corner(1)-bottom_right_corner(1)))/abs((bottom_left_corner(2)-bottom_right_corner(2))));
                toleft_angle = atand(abs((bottom_left_corner(1)-upper_left_corner(1)))/abs((bottom_left_corner(2)-upper_left_corner(2))));
                translate=[];

                if(toright_angle <= toleft_angle)
                    handles.rotation_angle = toright_angle;

                    % calculate the width and height of our rectangle
                    width = sqrt((bottom_left_corner(1)-bottom_right_corner(1))^2 ...
                        + (bottom_left_corner(2)-bottom_right_corner(2))^2);
                    height = sqrt((bottom_left_corner(1)-upper_left_corner(1))^2 ...
                        + (bottom_left_corner(2)-upper_left_corner(2))^2);

                    translate=bottom_left_corner;
                else
                    handles.rotation_angle = -toleft_angle;

                    % calculate the width and height of our rectangle
                    width = sqrt((bottom_left_corner(1)-upper_left_corner(1))^2 ...
                        + (bottom_left_corner(2)-upper_left_corner(2))^2);
                    height = sqrt((bottom_left_corner(1)-bottom_right_corner(1))^2 ...
                        + (bottom_left_corner(2)-bottom_right_corner(2))^2);

                    translate=upper_left_corner;
                end
                
                % save height and width
                handles.height = height;
                handles.width = width;
                
                % adjust sampling box determined by pixelspacing
                samplerate = get(handles.pixel_spacing,'String');  
                spacing = handles.info.PixelSpacing;
                spacing_factorx = spacing(2)/str2double(samplerate);
                spacing_factory = spacing(1)/str2double(samplerate);
                pseudo_height = height*spacing_factory;
                pseudo_width = width*spacing_factorx;
                
                % determine the number of pixels in our resulting ROI 
                % by rounding the width and length
                result = ones(round(pseudo_height),round(pseudo_width));
                %result2 = ones(round(height),round(width));
                %result3 = ones(round(height),round(width));

                %initialize rotation data array as emptys
                X=(1:size(handles.image,2));
                Y=(1:size(handles.image,1));
                Z=handles.image;

                % rotation matrix
                r=[cosd(handles.rotation_angle),sind(handles.rotation_angle);-sind(handles.rotation_angle),cosd(handles.rotation_angle)];     
                totalsize = size(result,1)*size(result,2);
                
                
                new_pointsx = result;
                new_pointsy = result;
                for i=1:size(result,1),
                    for j=1:size(result,2),
                        new_point = [1/spacing_factory*i;1/spacing_factorx*j];
                        new_point = r*new_point+translate;
                        new_pointsx(i,j) = new_point(2);
                        new_pointsy(i,j) = new_point(1);
                    end
                end
                
                result = interp2(X,Y,Z,new_pointsx,new_pointsy,'bicubic');
  
                handles.result = result;                       
            end
            handles.preview = figure;
            imshow(handles.result,'DisplayRange',[],'InitialMagnification',300);
        end
        handles.was_previewed = 1;
        set(handles.status_text,'String','Click confirm to save your ROI, or reset to start over.');
    else
        errordlg('You must first draw the ROI before you can preview!');
    end
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function ui_openfile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.just_filename,pathname] = uigetfile('*','Open a DICOM file',handles.opendir);
filename = strcat(pathname, handles.just_filename);
if isequal(filename,0) || isequal(pathname,0)
  errordlg('Please select a DICOM file');
else
  handles.opendir = pathname;
  handles.filename = filename;
  handles.image = double(dicomread(handles.filename));
  handles.info = dicominfo(handles.filename);
  handles.dim=size(handles.image);
  handles.cropped=handles.image; % let's make a copy
  handles.result=handles.image;
  handles.figure_image=imshow(handles.image,'DisplayRange',[]);
  handles.output = hObject;
  set(handles.status_text,'String','Click Draw ROI to select an ROI region, or click Load ROI to colocalize an ROI from another image.');
  handles.file_open = 1;
  if(~strcmp(handles.info.PatientName.FamilyName,handles.current_patient))
      handles.current_patient = handles.info.PatientName.FamilyName;
      set(handles.sequence_number,'String',num2str(1));
      set(handles.roi_number,'String',num2str(1));
      if(~handles.userset)
        handles.save_directory=cd;
      end
        save_text = strcat(handles.save_directory,'\',handles.info.PatientName.FamilyName);
        set(handles.save_directory_text,'String',save_text);
  end
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_openfile_Callback(hObject, eventdata, handles)
% hObject    handle to menu_openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.just_filename,pathname] = uigetfile('*','Open a DICOM file',handles.opendir);
filename = strcat(pathname, handles.just_filename);
if isequal(filename,0) || isequal(pathname,0)
  errordlg('Please select a DICOM file');
else
  handles.opendir = pathname;
  handles.filename = filename;
  handles.image = double(dicomread(handles.filename));
  handles.info = dicominfo(handles.filename);
  handles.dim=size(handles.image);
  handles.cropped=handles.image; % let's make a copy
  handles.result=handles.image;
  handles.figure_image=imshow(handles.image,'DisplayRange',[]);
  handles.output = hObject;
  set(handles.status_text,'String','Click Draw ROI to select an ROI region, or click Load ROI to colocalize an ROI from another image.');
  handles.file_open = 1;
  if(~strcmp(handles.info.PatientName.FamilyName,handles.current_patient))
      handles.current_patient = handles.info.PatientName.FamilyName;
      set(handles.sequence_number,'String',num2str(1));
      set(handles.roi_number,'String',num2str(1));
      if(~handles.userset)
        handles.save_directory=cd;
      end
        save_text = strcat(handles.save_directory,'\',handles.info.PatientName.FamilyName);
        set(handles.save_directory_text,'String',save_text);
  end
end
guidata(hObject,handles);

function status_text_Callback(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of status_text as text
%        str2double(get(hObject,'String')) returns contents of status_text as a double


% --- Executes during object creation, after setting all properties.
function status_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function contrast_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
 imcontrast(handles.figure_image);

 guidata(hObject,handles);
 
% --- Executes on button press in confirm.
function confirm_Callback(hObject, eventdata, handles)
% hObject    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw or load an ROI before attempting to confirm.');
else
    if(handles.roi_defined)
        handles.roi_defined =0;
        if(~handles.was_previewed)
            if(handles.is_rectangle)

                % save height and width
                handles.height = handles.y_max - handles.y_min;
                handles.width = handles.x_max - handles.x_min;

                % adjust sampling box determined by pixelspacing
                spacing = handles.info.PixelSpacing;
                spacing_factorx = spacing(2)/.5;
                spacing_factory = spacing(1)/.5;
                pseudo_height = handles.height*spacing_factory;
                pseudo_width = handles.width*spacing_factorx;

                % determine the number of pixels in our resulting ROI 
                % by rounding the width and length
                result = ones(round(pseudo_height),round(pseudo_width));
                translate = [handles.y_min;handles.x_min];

                %initialize rotation data array as emptys
                X=(1:size(handles.image,2));
                Y=(1:size(handles.image,1));
                Z=handles.image;         

                new_pointsx = result;
                new_pointsy = result;
                for i=1:size(result,1),
                    for j=1:size(result,2),
                        new_point = [1/spacing_factory*i;1/spacing_factorx*j];
                        new_point = new_point+translate;
                        new_pointsx(i,j) = new_point(2);
                        new_pointsy(i,j) = new_point(1);
                    end
                end

                result = interp2(X,Y,Z,new_pointsx,new_pointsy,'bicubic');

                handles.result = result;     
            else
                % variables used for rotation
                corners_x = handles.x_vals;
                corners_y = handles.y_vals;

                %find the bottom-left corner of the rectangle
                sorted_corners_y = sort(corners_y);
                sorted_corners_x = sort(corners_x);

                % find the lowest y coordinate
                lowest_corner = sorted_corners_y(1);        
                leftmost_corner = sorted_corners_x(4);
                bottom_left_corner = [lowest_corner;leftmost_corner];

                % find the lowest x coordinate that matches with the lowest y coordinate
                for i=1:4;                                  
                    if(corners_y(i) == lowest_corner && corners_x(i) <= leftmost_corner)
                        leftmost_corner = corners_x(i);
                        bottom_left_corner = [lowest_corner;leftmost_corner];
                    end
                end

                % find right-most lowest bottom corner
                rightmost_corner = sorted_corners_x(4);
                lowest_corner = sorted_corners_y(4);
                bottom_right_corner = [lowest_corner;rightmost_corner];
                for i=1:4;
                    if(corners_x(i) == rightmost_corner && corners_y(i) <= lowest_corner)
                        lowest_corner = corners_y(i);
                        bottom_right_corner = [lowest_corner;rightmost_corner];
                    end
                end

                % find left most corner not equal to bottom left corner
                leftmost_corner = sorted_corners_x(1);
                upper_left_corner = [0;leftmost_corner];
                for i=1:4;
                    if(corners_x(i) == leftmost_corner && corners_y(i) ~= bottom_left_corner(1))
                        upper_left_corner = [corners_y(i);leftmost_corner];
                    end
                end

                % set rotation angle
                toright_angle = atand(abs((bottom_left_corner(1)-bottom_right_corner(1)))/abs((bottom_left_corner(2)-bottom_right_corner(2))));
                toleft_angle = atand(abs((bottom_left_corner(1)-upper_left_corner(1)))/abs((bottom_left_corner(2)-upper_left_corner(2))));
                translate=[];

                if(toright_angle <= toleft_angle)
                    handles.rotation_angle = toright_angle;

                    % calculate the width and height of our rectangle
                    width = sqrt((bottom_left_corner(1)-bottom_right_corner(1))^2 ...
                        + (bottom_left_corner(2)-bottom_right_corner(2))^2);
                    height = sqrt((bottom_left_corner(1)-upper_left_corner(1))^2 ...
                        + (bottom_left_corner(2)-upper_left_corner(2))^2);

                    translate=bottom_left_corner;
                else
                    handles.rotation_angle = -toleft_angle;

                    % calculate the width and height of our rectangle
                    width = sqrt((bottom_left_corner(1)-upper_left_corner(1))^2 ...
                        + (bottom_left_corner(2)-upper_left_corner(2))^2);
                    height = sqrt((bottom_left_corner(1)-bottom_right_corner(1))^2 ...
                        + (bottom_left_corner(2)-bottom_right_corner(2))^2);

                    translate=upper_left_corner;
                end

                % save height and width
                handles.height = height;
                handles.width = width;

                % adjust sampling box determined by pixelspacing
                samplerate = get(handles.pixel_spacing,'String');  
                spacing = handles.info.PixelSpacing;
                spacing_factorx = spacing(2)/str2double(samplerate);
                spacing_factory = spacing(1)/str2double(samplerate);
                pseudo_height = height*spacing_factory;
                pseudo_width = width*spacing_factorx;

                % determine the number of pixels in our resulting ROI 
                % by rounding the width and length
                result = ones(round(pseudo_height),round(pseudo_width));

                %initialize rotation data array as emptys
                X=(1:size(handles.image,2));
                Y=(1:size(handles.image,1));
                Z=handles.image;

                % rotation matrix
                r=[cosd(handles.rotation_angle),sind(handles.rotation_angle);-sind(handles.rotation_angle),cosd(handles.rotation_angle)];     
                totalsize = size(result,1)*size(result,2);           

                new_pointsx = result;
                new_pointsy = result;
                for i=1:size(result,1),
                    for j=1:size(result,2),
                        new_point = [1/spacing_factory*i;1/spacing_factorx*j];
                        new_point = r*new_point+translate;
                        new_pointsx(i,j) = new_point(2);
                        new_pointsy(i,j) = new_point(1);
                    end
                end

                result = interp2(X,Y,Z,new_pointsx,new_pointsy,'bicubic');

                handles.result = result;     
            end
        end
        
        % variables used for saving
        original = handles.image;
        mask = handles.mask;
        metadata = handles.info;
        roi = double(handles.result);
        a_x = handles.a_x;
        a_y = handles.a_y;
        b_x = handles.b_x;
        b_y = handles.b_y;
        c_x = handles.c_x;
        c_y = handles.c_y;
        d_x = handles.d_x;
        d_y = handles.d_y;
        height = double(handles.height)*double(handles.info.PixelSpacing(1));
        width = double(handles.width)*double(handles.info.PixelSpacing(2));

        % change to the directory we want to save in
        save_directory = get(handles.save_directory_text,'String');            
        if(~exist(save_directory,'dir'))
            mkdir(save_directory);
        end      

        cd(save_directory);

        roi_number_str = get(handles.roi_number,'String');
        sequence_number_str = get(handles.sequence_number,'String');
        roi_number_num = str2num(roi_number_str);
        sequence_number_num = str2num(sequence_number_str);
        patient_name = handles.info.PatientName.FamilyName;
        if(strcmp(patient_name,''))
            patient_name = 'UNKNOWN_PATIENT';
        end

        output_filename = strcat(patient_name,'_SEQUENCE',sequence_number_str,'_ROI',roi_number_str,'.mat');

        save(output_filename, 'roi','mask','original','metadata','a_x','a_y','b_x','b_y','c_x','c_y','d_x','d_y','height','width'); 
        output_dialog = sprintf('Your ROI has been saved as %s', output_filename);
        set(handles.status_text,'String',output_dialog);

        % reset
            
        v = caxis;
        r = axis;

        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.was_previewed=0;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);

        % update the sequence number
        if(roi_number_num == 5 && sequence_number_num ~= 5)
            set(handles.sequence_number,'String',num2str(sequence_number_num + 1));
        elseif(roi_number_num == 5 && sequence_number_num == 5)
            set(handles.sequence_number,'String',num2str(1));
        end

        % update the roi number
        if(roi_number_num == 5)
            set(handles.roi_number,'String',num2str(1));
        else
            set(handles.roi_number,'String',num2str(roi_number_num + 1));
        end

        % change back to the workspace directory
        cd(handles.workspace);

        guidata(hObject,handles);
    else
        errordlg('You must first draw an ROI before attempting to confirm!');
    end

end


% --- Executes on button press in load_roi.
function load_roi_Callback(hObject, eventdata, handles)
% hObject    handle to load_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file before attempting to load an ROI.');
else
    if(handles.roi_defined)    
        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        handles.was_previewed=0;
    end
    
    matching_patient = 0;
    
    while(~matching_patient)        
        set(handles.status_text,'String','Select a file to define a colocalized ROI.');
        [filename,pathname] = uigetfile('*','Open a MAT file',handles.loaddir);
        handles.loaddir=pathname;
        handles.mat_file = strcat(pathname,filename);
        if isequal(filename,0) || isequal(pathname,0)
          errordlg('Please select a MAT file')
        else
          load(handles.mat_file);
          handles.original = original;
          previousinfo = metadata;
          handles.use_load = 1;
          handles.roi = roi;
          handles.mask = mask;
          handles.a_x = a_x;
          handles.a_y = a_y;
          handles.b_x = b_x;
          handles.b_y = b_y;
          handles.c_x = c_x;
          handles.c_y = c_y;
          handles.d_x = d_x;
          handles.d_y = d_y;
          set(handles.status_text,'String','Move the ROI as desired before clicking Preview.');
          
          if previousinfo.PatientName.FamilyName == handles.info.PatientName.FamilyName
              matching_patient=1;
          end
        end
    
        if(~matching_patient)
            e = errordlg('You must choose a file which matches the current patient. Try again.');
            waitfor(e);
        end
    end

    m_atob = (handles.a_y - handles.b_y)/(handles.a_x - handles.b_x);
    m_btoc = -1/m_atob;

    atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
    handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;
    btoc_theta_inradians = atan(m_btoc);
    handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;

    % special case, handle perpendicular rectangles
    if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
        % we have a normal non-rotated rectangle
        handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
        handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
        handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
        handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

        rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
        handles.is_rectangle = 1;

        handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
        handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
        handles.mask = roipoly(handles.image, handles.X, handles.Y);


    else
        handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
        handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

        handles.line = line(handles.X,handles.Y,'Marker','.','LineStyle','-');
        handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
        handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];

        handles.mask = roipoly(handles.image, handles.X, handles.Y);
    end

    handles.roi_defined =1;
    set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
    % Update handles structure
    guidata(hObject,handles);
end


% --- Executes on button press in move_up.
function move_up_Callback(hObject, eventdata, handles)
% hObject    handle to move_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw an ROI before attempting to move it.');
else
    if(~handles.roi_defined)
        errordlg('You must first define an ROI before you can move it.')
    else
        % get contrast/zoom settings
        v = caxis;
        r = axis;
        
        % reset image
        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        handles.was_previewed=0;
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);
 

        handles.a_y = handles.a_y - .5;
        handles.b_y = handles.b_y - .5;
        handles.c_y = handles.c_y - .5;
        handles.d_y = handles.d_y - .5;

        % find slopes
        m_atob = (handles.a_y - handles.b_y)/(handles.a_x - handles.b_x);
        m_btoc = -1/m_atob;
        m_ctod = m_atob;
        m_dtoa = m_btoc;        

        atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
        handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;
        btoc_theta_inradians = atan(m_btoc);
        handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;
       

        % special case, handle perpendicular rectangles
        if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
            % we have a normal non-rotated rectangle
            handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
            handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

            rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
            handles.is_rectangle = 1;

            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
            handles.mask = roipoly(handles.image, handles.X, handles.Y);


        else
            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

            line(handles.X,handles.Y,'Marker','.','LineStyle','-');
            handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
            handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];
            handles.mask = roipoly(handles.image, handles.x_vals, handles.y_vals);

            handles.imagecopy = imcrop(handles.image,[min(handles.x_vals)-(max(handles.x_vals)-min(handles.x_vals)),min(handles.y_vals)-(max(handles.y_vals)-min(handles.y_vals)),3*(max(handles.x_vals)-min(handles.x_vals)),3*(max(handles.y_vals)-min(handles.y_vals))]);
            handles.ones = ones(size(handles.image));
            handles.ones(round(handles.a_y),round(handles.a_x)) = 255;
            handles.ones(round(handles.b_y),round(handles.b_x)) = 255;
            handles.ones(round(handles.c_y),round(handles.c_x)) = 255;
            handles.ones(round(handles.d_y),round(handles.d_x)) = 255;
            handles.ones = imcrop(handles.ones,[min(handles.x_vals)-(max(handles.x_vals)-min(handles.x_vals)),min(handles.y_vals)-(max(handles.y_vals)-min(handles.y_vals)),3*(max(handles.x_vals)-min(handles.x_vals)),3*(max(handles.y_vals)-min(handles.y_vals))]);

        end


        handles.roi_defined =1;
        set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
        % Update handles structure
        guidata(hObject, handles);
    end
end
    

% --- Executes on button press in move_down.
function move_down_Callback(hObject, eventdata, handles)
% hObject    handle to move_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw an ROI before attempting to move it.');
else
    if(~handles.roi_defined)
        errordlg('You must first define an ROI before you can move it.')
    else
        % get contrast/zoom settings
        v = caxis;
        r = axis;
        
        % reset image
        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        handles.was_previewed=0;
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);

        handles.a_y = handles.a_y + .5;
        handles.b_y = handles.b_y + .5;
        handles.c_y = handles.c_y + .5;
        handles.d_y = handles.d_y + .5;

        % calculate slopes
        m_atob = (handles.a_y - handles.b_y)/(handles.a_x - handles.b_x);
        m_btoc = -1/m_atob;
        m_ctod = m_atob;
        m_dtoa = m_btoc;


        atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
        handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;
        btoc_theta_inradians = atan(m_btoc);
        handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;

        % special case, handle perpendicular rectangles
        if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
            % we have a normal non-rotated rectangle
            handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
            handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

            rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
            handles.is_rectangle = 1;

            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
            handles.mask = roipoly(handles.image, handles.X, handles.Y);


        else
            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

            line(handles.X,handles.Y,'Marker','.','LineStyle','-');
            handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
            handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];
            handles.mask = roipoly(handles.image, handles.x_vals, handles.y_vals);

        end


        handles.roi_defined =1;
        set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --- Executes on button press in move_left.
function move_left_Callback(hObject, eventdata, handles)
% hObject    handle to move_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw an ROI before attempting to move it.');
else
    if(~handles.roi_defined)
        errordlg('You must first define an ROI before you can move it.');
    else
        % get contrast/zoom settings
        v = caxis;
        r = axis;
        
        % reset image
        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        handles.was_previewed=0;
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);
        
        handles.a_x = handles.a_x - .5;
        handles.b_x = handles.b_x - .5;
        handles.c_x = handles.c_x - .5;
        handles.d_x = handles.d_x - .5;

        % slopes
        m_atob = (handles.a_y - handles.b_y)/(handles.a_x - handles.b_x);
        m_btoc = -1/m_atob;
        m_ctod = m_atob;
        m_dtoa = m_btoc;

        atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
        handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;
        btoc_theta_inradians = atan(m_btoc);
        handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;

        % special case, handle perpendicular rectangles
        if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
            % we have a normal non-rotated rectangle
            handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
            handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

            rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
            handles.is_rectangle = 1;

            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
            handles.mask = roipoly(handles.image, handles.X, handles.Y);


        else
            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

            line(handles.X,handles.Y,'Marker','.','LineStyle','-');
            handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
            handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];
            handles.mask = roipoly(handles.image, handles.x_vals, handles.y_vals);

        end


        handles.roi_defined =1;
        set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --- Executes on button press in move_right.
function move_right_Callback(hObject, eventdata, handles)
% hObject    handle to move_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.file_open)
    errordlg('Please open a file and draw an ROI before attempting to move it.');
else
    if(~handles.roi_defined)
        errordlg('You must first define an ROI before you can move it.')
    else
        % get contrast/zoom settings
        v = caxis;
        r = axis;
        
        % reset image
        handles.image = dicomread(handles.filename);
        handles.dim=size(handles.image);
        handles.cropped=handles.image; % let's make a copy
        handles.result=handles.image;
        handles.figure_image=imshow(handles.image,'DisplayRange',[]);
        handles.output = hObject;
        handles.was_previewed=0;
        caxis([v(1) v(2)]);
        axis([r(1) r(2) r(3) r(4)]);
        
        handles.a_x = handles.a_x + .5;
        handles.b_x = handles.b_x + .5;
        handles.c_x = handles.c_x + .5;
        handles.d_x = handles.d_x + .5;

        % slopes
        m_atob = (handles.a_y - handles.b_y)/(handles.a_x - handles.b_x);
        m_btoc = -1/m_atob;
        m_ctod = m_atob;
        m_dtoa = m_btoc;

        atob_theta_inradians = atan(sqrt((handles.a_y-handles.b_y)^2)/sqrt((handles.a_x-handles.b_x)^2));
        handles.atob_theta_indegrees = 180 * atob_theta_inradians / pi;
        btoc_theta_inradians = atan(m_btoc);
        handles.btoc_theta_indegrees = 180 * btoc_theta_inradians / pi;

        % calculate y-intercepts
        b_atob = handles.b_y - m_atob*handles.b_x;
        b_btoc = handles.c_y - m_btoc*handles.c_x;
        b_ctod = handles.d_y - m_ctod*handles.d_x;
        b_dtoa = handles.a_y - m_dtoa*handles.a_x;

        % special case, handle perpendicular rectangles
        if(handles.a_x == handles.b_x || handles.a_y == handles.b_y || handles.b_x == handles.c_x || handles.b_y == handles.c_y)
            % we have a normal non-rotated rectangle
            handles.x_min = min([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_min = min([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);
            handles.x_max = max([handles.a_x,handles.b_x,handles.c_x,handles.d_x]);
            handles.y_max = max([handles.a_y,handles.b_y,handles.c_y,handles.d_y]);

            rectangle('Position',[handles.x_min,handles.y_min,handles.x_max-handles.x_min,handles.y_max-handles.y_min],'EdgeColor', 'blue');
            handles.is_rectangle = 1;

            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];
            handles.mask = roipoly(handles.image, handles.X, handles.Y);

        else
            handles.X = [handles.a_x handles.b_x handles.c_x handles.d_x handles.a_x];
            handles.Y = [handles.a_y handles.b_y handles.c_y handles.d_y handles.a_y];

            line(handles.X,handles.Y,'Marker','.','LineStyle','-');
            handles.x_vals = [handles.a_x, handles.b_x, handles.c_x, handles.d_x];
            handles.y_vals = [handles.a_y, handles.b_y, handles.c_y, handles.d_y];
            handles.mask = roipoly(handles.image, handles.x_vals, handles.y_vals);

        end

        handles.roi_defined =1;
        set(handles.status_text,'String','Hit Reset to redraw your ROI, Preview to view it, or Confirm to crop the ROI.');
        % Update handles structure
        guidata(hObject, handles);
    end
end



function sequence_number_Callback(hObject, eventdata, handles)
% hObject    handle to sequence_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sequence_number as text
%        str2double(get(hObject,'String')) returns contents of sequence_number as a double


% --- Executes during object creation, after setting all properties.
function sequence_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sequence_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roi_number_Callback(hObject, eventdata, handles)
% hObject    handle to roi_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_number as text
%        str2double(get(hObject,'String')) returns contents of roi_number as a double


% --- Executes during object creation, after setting all properties.
function roi_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;


% --- Executes on button press in save_directory.
function save_directory_Callback(hObject, eventdata, handles)
% hObject    handle to save_directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.save_directory= uigetdir('*','Choose the base directory you would like to save your files in...'); 
if isequal(handles.save_directory,0)
  errordlg('Please select a directory.')
else    
    if(~handles.file_open)
        save_text = strcat(handles.save_directory,'\');
    else
        save_text = strcat(handles.save_directory,'\',handles.info.PatientName.FamilyName); 
    end
    set(handles.save_directory_text,'String',save_text);
    handles.userset=1;
end
guidata(hObject,handles);



function save_directory_text_Callback(hObject, eventdata, handles)
% hObject    handle to save_directory_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_directory_text as text
%        str2double(get(hObject,'String')) returns contents of save_directory_text as a double


% --- Executes during object creation, after setting all properties.
function save_directory_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_directory_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
% hObject    handle to Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function howto_Callback(hObject, eventdata, handles)
% hObject    handle to howto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = strvcat('ROIDraw GUI Help',...
' ',...
'File: Open - opens a new DICOM file.',...
'File: Quit - quits the program.',...
' ',...
'ICONS:',...
'Folder Icon - opens a new DICOM file.',...
'Magnifying Glass(+) Icon - increases magnification at a specified location.',...
'Magnifying Glass(-) Icon- decreases magnification at a specified location.',...
'Cursor Icon - gives information of pixel at a specified location.',...
'Hand Icon - allows you to drag the image around to view different parts of the image.',...
'Black/White Circle Icon - allows you to alter the contrast/brightness settings.',...
' ',...
'STATUS DIALOGS',...
'Left Dialog - gives you tips on how to use ROIDraw.',...
'Right Dialog - shows the directory in which your ROI will be saved. You can edit this by hand if you would like.',...
' ',...
'DRAWING CONTROLS',...
'Arrows - allows you to move a drawn ROI up/down and left/right.',...
'Load ROI - opens a dialog to allow you to open a file from a previously drawn and saved ROI and places an ROI at the same location.',...
'Reset - resets the screen. Manigification and contrast settings are kept however.',...
'Preview - opens a new window to view the ROI you have selected.',...
' ',...
'SAVING CONTROLS',...
'Sample Rate: allows you to control the rate at which ROI's are sampled from the original image.',...
'Sequence Number: shows the sequence number that will be in the saved file name. You can edit this yourself.',...
'ROI Number: shows the ROI number that will be in the saved file name. You can edit this yourself. ',...
'Save Directory: gives you control over which base directory the files will be saved in. Files will be saved in a subfolder by patient name unless you alter the directory by hand in the dialog box.',...
'Confirm: saves your ROI and resets the image',...
'*Note: Files are all saved as PATIENTNAME_SEQUENCE#_ROI#.mat'...
);
HelpBox = helpdlg(str,'HELP');



function pixel_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_spacing as text
%        str2double(get(hObject,'String')) returns contents of pixel_spacing as a double


% --- Executes during object creation, after setting all properties.
function pixel_spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


