IM_seq = {'./data/MainFrame_1.png',...
'./data/MainFrame_2.png',...
'./data/MainFrame_3.png',...
'./data/MainFrame_4.png',...
'./data/MainFrame_5.png',...
'./data/MainFrame_6.png',...
'./data/MainFrame_7.png',...
'./data/MainFrame_8.png',...
'./data/MainFrame_9.png',...
'./data/MainFrame_10.png',...
'./data/MainFrame_11.png',...
'./data/MainFrame_12.png',...
'./data/MainFrame_13.png',...
'./data/MainFrame_14.png',...
'./data/MainFrame_23.png',...
'./data/MainFrame_24.png',...
'./data/MainFrame_25.png',...
'./data/MainFrame_26.png',...
'./data/MainFrame_27.png',...
'./data/MainFrame_28.png',...
'./data/MainFrame_29.png'};
SE_seq = {'./data/Segmentation_1.png',...
'./data/Segmentation_2.png',...
'./data/Segmentation_3.png',...
'./data/Segmentation_4.png',...
'./data/Segmentation_5.png',...
'./data/Segmentation_6.png',...
'./data/Segmentation_7.png',...
'./data/Segmentation_8.png',...
'./data/Segmentation_9.png',...
'./data/Segmentation_10.png',...
'./data/Segmentation_11.png',...
'./data/Segmentation_12.png',...
'./data/Segmentation_13.png',...
'./data/Segmentation_14.png',...
'./data/Segmentation_23.png',...
'./data/Segmentation_24.png',...
'./data/Segmentation_25.png',...
'./data/Segmentation_26.png',...
'./data/Segmentation_27.png',...
'./data/Segmentation_28.png',...
'./data/Segmentation_29.png'};

% load image pakge for image operations
pkg load image

IM = {};
SE = {};
for iIM = 1:numel(IM_seq)
  IM{iIM} = imread(IM_seq{iIM});
  SE{iIM} = imread(SE_seq{iIM});
end

IM1 = IM{1};
IM2 = SE{1};

IMblend = zeros(size(IM1));

segmentationcolors = [255 0 0; 0 255 0; 0 0 255; 0 0 0; 255 255 0; 0 255 255; 255 255 255; 128 128 128];
segmentationnames = {'Hair'; 'Skin'; 'Eyes'; 'Torso'; 'Arms'; 'Belt'; 'Pants'; 'Shoes'};
h.segmentationnames = segmentationnames;
h.segmentationcolors = segmentationcolors;

sliderdefault = [.22, .14, .64, .1, .5, .5, .66, .5];
blendingdefault = [3, 2, 2, 2, 3, 3, 2, 3];
blendingmodes = {'g'; 'o'; 'o'; 'o'; 'g'; 'g'; 'o'; 'g'}; % m = multiply; o=overlay; g=gray
h.bmodes = {'m','o','g'};
h.imgs = IM;
h.segs = SE;
h.cutrows = [];
h.cutcols = [];
h.background = [255  0 255];

comap_rainbow = [[(0:0.1:1)', [0:0.033:0.16, 0.18:-0.033:0]', zeros(11,1)];rainbow(234);[(2/3:1/30:1)', (0:0.1:1)', ones(11,1)]]; % custom rainbow colormap with a bit of black and white on the edges
comap_hsv = [[(0:0.1:1)', [0:0.033:0.16, 0.18:-0.033:0]', zeros(11,1)];hsv(234); [ones(11,1), (0:0.1:1)', (0.0256:0.09:1)']];

%% explanation of blending modes:
% multiply: I1*I2 / 255
% overlay: I1/255 * (I1 + (2*I2)*(255-I1)/255)

%%
%% hair = red
%% skin = green
%% eyes = blue
%% torso = black
%% arms = yellow
%% belt = cyan
%% pants = white
%% shoe = grey

%% ---- callback functions
function load_sequence(obj, init = false)
  % open file dialog for frames
  [fname_frame, fpath_frame, fltidx_frame] = uigetfile ({"*.png", "Frames in png format"},"Load frames","*.png", "Multiselect", "on");
  h = guidata (obj);
  bConv = logical(get(h.BCb, "value"));
  % read frames
  if ~isnumeric(fname_frame)
    if iscell(fname_frame)
      for iIM = 1:numel(fname_frame)
        % convert rgb2gray if required
        if bConv
          IM_tmp = imread(fullfile(fpath_frame,fname_frame{iIM}));
          IM_tmp_bg = (IM_tmp(:,:,1) == 255) & (IM_tmp(:,:,2) == 0) & (IM_tmp(:,:,3) == 255); % get magenta background
          IM_tmp_ = repmat(rgb2gray(IM_tmp),1,1,3);
          IM_tmp(~repmat(IM_tmp_bg,1,1,3)) = IM_tmp_(~repmat(IM_tmp_bg,1,1,3));
          IM{iIM} = IM_tmp;
        else
          IM{iIM} = imread(fullfile(fpath_frame,fname_frame{iIM}));
        end
      end
    else
      % convert rgb2gray if required
      if bConv
        IM_tmp = imread(fullfile(fpath_frame,fname_frame));
        IM_tmp_bg = (IM_tmp(:,:,1) == 255) & (IM_tmp(:,:,2) == 0) & (IM_tmp(:,:,3) == 255); % get magenta background
        IM_tmp_ = repmat(rgb2gray(IM_tmp),1,1,3);
        IM_tmp(~repmat(IM_tmp_bg,1,1,3)) = IM_tmp_(~repmat(IM_tmp_bg,1,1,3));
        IM{1} = IM_tmp;
      else
        IM{1} = imread(fullfile(fpath_frame,fname_frame));
      end
    end
  else
    warning("No frames selected.")
    return;
  end
  % open file dialog for segmentation maps
  [fname_seg, fpath_seg, fltidx_seg] = uigetfile ({"*.png", "Segmentations in png format"},"Load Segmentation Maps","*.png", "Multiselect", "on");
  % read segmentation maps
  if ~isnumeric(fname_seg)
    if iscell(fname_seg)
      for iSE = 1:numel(fname_seg)
        SE{iSE} = imread(fullfile(fpath_seg,fname_seg{iSE}));
      end
    else
      SE{1} = imread(fullfile(fpath_seg,fname_seg));
    end
  else
    warning("No segmentation maps selected.")
    return;
  end
  % check that number of images is equal to number of frame and all have equal size
  if numel(IM) ~= numel(SE)
    warning("Number of frames and segmentation maps must be equal.");
    return;
  end
  % update guidata
  h.imgs = IM;
  h.segs = SE;
  set(h.PS, "string", arrayfun(@num2str,1:numel(IM),'UniformOutput',false));
  set(h.PS, "value", 1);
  guidata(obj, h);
  % update cut and image
  update_cut(obj, true);
end

% check if rows or colums should be cut from input
function update_cut(obj, init = false)
  h = guidata (obj);
  PCin = str2num(get(h.PC, "string"));
  PCcin = str2num(get(h.PCc, "string"));
  if ~isempty(PCin)
    if isnumeric(PCin)
      h.cutrows = unique(PCin);
    end
  else
    h.cutrows = [];
  end
  if ~isempty(PCcin)
    if isnumeric(PCcin)
      h.cutcols = unique(PCcin);
    end
  else
    h.cutcols = [];
  end
  guidata(obj, h);
  update_image(obj, true);
end

% perform image operations and update view
function update_image(obj, init = false)
  h = guidata (obj);
  sIDX = get(h.PS, "value");
  dVal = (get(h.LGE, "value")-1);
  % get original segmentation maps
  masks_orig = CreateMasks(h.segs{sIDX}, h.segmentationcolors);
  % grow hair (Experimental)
  IM_tmp = h.imgs{sIDX};
  h.IMin = h.imgs{sIDX};
  SEGs = h.segs{sIDX};
  % dilate hair mask
  if dVal > 0
    SE = strel("disk",ceil(dVal),0);
    % do stuff
    IM_tmp(:,:,1) = 0;
    IM_tmp(:,:,2) = 0;
    IM_tmp(:,:,3) = 0;
    IM_tmp(masks_orig{1}) = h.IMin(masks_orig{1});
    IM_tmp = imdilate(IM_tmp, SE);
    % grow hair mask
    hair_mask_dlt = imdilate(masks_orig{1},SE);
    h.IMin(hair_mask_dlt) = IM_tmp(hair_mask_dlt);
    % modify segmentations (red) to grow along hair
    SEGs_red = SEGs(:,:,1);
    SEGs_green = SEGs(:,:,2);
    SEGs_blue = SEGs(:,:,3);
    SEGs_red(hair_mask_dlt(:,:,1)) = 255;
    SEGs_green(hair_mask_dlt(:,:,2)) = 0;
    SEGs_blue(hair_mask_dlt(:,:,3)) = 0;
    SEGs(:,:,1) = SEGs_red;
    SEGs(:,:,2) = SEGs_green;
    SEGs(:,:,3) = SEGs_blue;
  end
  %
  % cut rows and columns
  if ~isempty(h.cutrows)
    h.IMin = cutrows(h.IMin,h.cutrows,h.background);
    masks_new = CreateMasks(cutrows(SEGs,h.cutrows,h.background), h.segmentationcolors);
  else
    %h.IMin = h.imgs{sIDX};
    masks_new = CreateMasks(SEGs, h.segmentationcolors);
  end
  if ~isempty(h.cutcols)
    h.IMin = cutcols(h.IMin,h.cutcols,h.background);
    masks_new = CreateMasks(cutcols(cutrows(SEGs,h.cutrows,h.background), h.cutcols, h.background), h.segmentationcolors);
  end
  h.masks = masks_new; %update masks
  guidata(obj, h); % update guidata
  update_plot(obj, true);   
end

% update view
function update_plot (obj, init = false)
   h = guidata (obj);
   blendingmodes = h.blendingmodes;
   blendingcolors = h.blendingcolors;
   % set all blending modes
   for iMode = 1:numel(h.blendingmodes)
     blendingmodes{iMode} = h.bmodes{get (h.(['p',num2str(iMode)]), "value")};
   end
   % set all blending colors
   for iMode = 1:numel(h.blendingmodes)
     svalue = get(h.(['s',num2str(iMode)]), "value");
     if blendingmodes{iMode} == 'g'
       blendingcolors(iMode,:) = [svalue svalue svalue];  
     else
       blendingcolors(iMode,:) = h.comap(round(svalue*255+1),:);
     end
     set(h.(['l',num2str(iMode)]), "string", num2str(svalue,2));
   end
   % blend image
   IMblend = blending(h.IMin, h.masks, blendingcolors, blendingmodes);
   set(h.IM, 'CData', IMblend);
end

% export current frame (CData) to png file
function export_current_frame (obj, init = false)
  h = guidata (obj);
  CData = get(h.IM, 'CData');
  [fname, fpath, fltidx] = uiputfile ({"*.png", "png format"},"Save frame to","Frame.png");
  if ~isnumeric(fname) & ~isnumeric(fpath)
    imwrite(CData, fullfile(fpath,fname));
  end
end

% export all frames to png files
function export_all_frames(obj, init = false)
  h = guidata (obj);
  current_frame = get(h.PS, "value");
  for iIM = 1:numel(h.imgs)
    set(h.PS, "value", iIM);
    update_cut(obj, true);
    drawnow;
    h = guidata (obj);
    CData = get(h.IM, 'CData');
    imwrite(CData, [get(h.EED, "string"), num2str(iIM), '.png']);
  end
  set(h.PS, "value", current_frame);
  h = guidata (obj);
  update_cut(obj, true);
end

% perform image operations blending
function IMout = blending(IMin, masks, blendingcolors, blendingmodes)
  UImat = ones(size(IMin));
  IMblend = IMin;
  for iPart = 1:size(blendingcolors,1)
    mask_tmp = masks{iPart};
    IM_red = double(IMin(:,:,1));
    IM_green = double(IMin(:,:,2));
    IM_blue = double(IMin(:,:,3));
    IMblend_red = double(IMblend(:,:,1));
    IMblend_green = double(IMblend(:,:,2));
    IMblend_blue = double(IMblend(:,:,3));
    % blending
    if blendingmodes{iPart} == 'm' %multiply I1*I2 / 255
      IMblend_red(mask_tmp(:,:,1)) = IM_red(mask_tmp(:,:,1))*blendingcolors(iPart,1);
      IMblend_green(mask_tmp(:,:,2)) = IM_green(mask_tmp(:,:,2))*blendingcolors(iPart,2);
      IMblend_blue(mask_tmp(:,:,3)) = IM_blue(mask_tmp(:,:,3))*blendingcolors(iPart,3);
    elseif blendingmodes{iPart} == 'o' %overlay I1/255 * (I1 + (2*I2)*(255-I1)/255)
      IMblend_red(mask_tmp(:,:,1)) = IM_red(mask_tmp(:,:,1))/255 .* (IM_red(mask_tmp(:,:,1)) + 2*UImat(mask_tmp(:,:,1))*blendingcolors(iPart,1) .* (255 - IM_red(mask_tmp(:,:,1))) );
      IMblend_green(mask_tmp(:,:,2)) = IM_green(mask_tmp(:,:,2))/255 .* (IM_green(mask_tmp(:,:,2)) + 2*UImat(mask_tmp(:,:,2))*blendingcolors(iPart,2) .* (255 - IM_green(mask_tmp(:,:,2))) );
      IMblend_blue(mask_tmp(:,:,3)) = IM_blue(mask_tmp(:,:,3))/255 .* (IM_blue(mask_tmp(:,:,3)) + 2*UImat(mask_tmp(:,:,3))*blendingcolors(iPart,3) .* (255 - IM_blue(mask_tmp(:,:,3))) );
    elseif blendingmodes{iPart} == 'g' %gray-tone gamma filter, .* (I/128)^1.5
      IMblend_red(mask_tmp(:,:,1)) = IM_red(mask_tmp(:,:,1))*(blendingcolors(iPart,1)*2)^1.5;
      IMblend_green(mask_tmp(:,:,2)) = IM_green(mask_tmp(:,:,2))*(blendingcolors(iPart,1)*2)^1.5;
      IMblend_blue(mask_tmp(:,:,3)) = IM_blue(mask_tmp(:,:,3))*(blendingcolors(iPart,1)*2)^1.5;
    end
    % recompose image
    IMblend(:,:,1) = uint8(IMblend_red);
    IMblend(:,:,2) = uint8(IMblend_green);
    IMblend(:,:,3) = uint8(IMblend_blue);
  end
IMout = IMblend;
end

% create segmentation marks from pre-defined segmentation colors
function masks_out = CreateMasks(IMseg, segmentationcolors)
  masks_out = {};
  for iPart = 1:size(segmentationcolors,1)
	  mask_tmp = (IMseg(:,:,1) == segmentationcolors(iPart,1)) & (IMseg(:,:,2) == segmentationcolors(iPart,2)) & (IMseg(:,:,3) == segmentationcolors(iPart,3));
	  masks_out{end+1} = repmat(mask_tmp,[1 1 3]);
  end
end

% cut rows from input image
function IMout = cutrows(IMin, rows, bgcolor)
  % create blank img
  IMout = ones(size(IMin), 'uint8');
  IMout(:,:,1) = bgcolor(1);
  IMout(:,:,2) = bgcolor(2);
  IMout(:,:,3) = bgcolor(3);
  % Cut rows
  IMin(rows,:,:) = [];
  % Shift rest of image down
  IMout((numel(rows)+1):end,:,:) = IMin;
end

% cut columns from input image
function IMout = cutcols(IMin, cols, bgcolor)
  % create blank img
  IMout = ones(size(IMin), 'uint8');
  IMout(:,:,1) = bgcolor(1);
  IMout(:,:,2) = bgcolor(2);
  IMout(:,:,3) = bgcolor(3);
  % Cut rows
  IMin(:,cols,:) = [];
  % Shift rest of image down
  IMout(:,(numel(cols)+ceil(numel(cols)/2)):end,:) = IMin;
end

%% --- do some checks maybe (are segmentation colors equal to colors used in IM2?)
if ~all(size(IM1) == size(IM2))
	error('Frame and segmentation mask need to have equal size.');
end

if ~(numel(IM_seq) == numel(SE_seq))
	error('Number of frames and segmentation maps must be identical.');
end

if ~all(size(segmentationcolors,1) == size(segmentationnames,1))
	error('Incorrect number of segmentations. Check number of segmentation colors and names.');
end

if ~all(size(segmentationcolors,1) == size(sliderdefault,2))
  error('Incorrect number of segmentations or default values.');
end

% Create 3-channel binary masks for each part
masks = {};
masks = CreateMasks(IM2, segmentationcolors);

h.IMin =  IM1;
h.masks = masks;
h.blendingmodes = blendingmodes;
h.comap = comap_hsv;
h.blendingcolors = [arrayfun(@(x) h.comap(round(x*255+1),1),sliderdefault')...
 arrayfun(@(x) h.comap(round(x*255+1),2),sliderdefault')...
 arrayfun(@(x) h.comap(round(x*255+1),3),sliderdefault')];

IMblend = blending(IM1, masks, h.blendingcolors, h.blendingmodes);

% draw the GUI
hf = figure('Position',[300 200 900 600]);
axIM = axes('units', 'normalized', 'Position',[0.1 0.0 .4 1.0]);
h.IM = image(IMblend);
axis image;
axis off;
p = uipanel ("title", "Colors", "position", [.5 .01 .4 .97]);

% plot a jet colorbar as reference
ax_co = axes( p, 'units', 'normalized', 'position', [0.05, 0.02, 0.9, 0.1] );
C = [];
C(:,:,2) = h.comap;
C = permute(C,[3 1 2]);
% add gray bar as top column
C(1,:,1) = (0:255)/255;
C(1,:,2) = (0:255)/255;
C(1,:,3) = (0:255)/255;
imagesc(ax_co, C);
axis off;

% add control elements to panel:
% 'hair'; 'skin'; 'eyes'; 'torso'; 'arms'; 'belt'; 'pants'; 'shoes'
% ToDo: Load sequence (frames and segmentation maps)
% ToDo: Export series as png
pos_top = .93;
pos_step = .108;
for iMasks = 1:numel(segmentationnames)
  h.(['t',num2str(iMasks)]) = uicontrol ("parent", p, "style", "text", "string", segmentationnames{iMasks}, 'units', 'normalized', "position",[.01 pos_top-(iMasks-1)*pos_step .15 .05]);
  h.(['l',num2str(iMasks)]) = uicontrol ("parent", p, "style", "text", "string", num2str(sliderdefault(iMasks)), 'units', 'normalized', "position",[.87 pos_top-(iMasks-1)*pos_step .15 .05]);
  h.(['p',num2str(iMasks)]) = uicontrol ("parent", p, "style", "popupmenu", "string", {"Multiply"; "Overlay"; "Gray"}, 'value', blendingdefault(iMasks), 'units', 'normalized',  "position", [.2 pos_top-(iMasks-1)*pos_step .3 .05], "callback", @update_plot);
  h.(['s',num2str(iMasks)]) = uicontrol ("parent", p, "style", "slider", "value", sliderdefault(iMasks) , 'units', 'normalized', "position",[.02 pos_top-(iMasks-1)*pos_step-.05 .95 .04], "callback", @update_plot);
end
% Grow edit for hair
h.LGR = uicontrol("parent", p, "style", "text", "string", "Grow:", 'units', 'normalized', "position",[.55 pos_top+0.01 .15 .04]);
h.LGE = uicontrol("parent", p, "style", "popupmenu", "string", {"0";"1";"2";"3"}, 'value', 1, 'units', 'normalized', "position",[.7 pos_top+0.01 .12 .04], "callback", @update_image);

% add open sequence button
h.BO = uicontrol("style", "pushbutton", "string", "Open frames", 'units', 'normalized', "position",[.005 .95 .12 .04], "callback", @load_sequence); 
h.BCb = uicontrol("style", "checkbox", "string", "rgb2gray", "value", 0, "units", "normalized", "position", [.005 .90 .12 .04]);

% add selection of frames
h.LF = uicontrol("style", "text", "string", "Frame:", 'units', 'normalized', "position",[.003 .80 .08 .04]); 
h.PS = uicontrol("style", "popupmenu", "string", arrayfun(@num2str,1:numel(IM),'UniformOutput',false), 'value', 1, 'units', 'normalized', "position",[.005 .75 .10 .05], "callback", @update_image); 

% Export frame
h.EBT = uicontrol("style", "pushbutton", "string", {'Save frame'}, 'units', 'normalized', "position",[.005 .65 .12 .06], "callback", @export_current_frame);
h.ESBT = uicontrol("style", "pushbutton", "string", {'Save all'}, 'units', 'normalized', "position",[.005 .57 .12 .06], "callback", @export_all_frames);
h.EED = uicontrol("style", "edit", "string", "Frame_edit_", 'units', 'normalized', "position",[.005 .53 .12 .04]);

% Cut row option
h.LC = uicontrol("style", "text", "string", "Cut rows:", 'units', 'normalized', "position",[.003 .75-0.3 .12 .04]); 
h.PC = uicontrol("style", "edit", "string", "", 'units', 'normalized', "position",[.005 .70-0.3 .12 .05]); 
h.BC = uicontrol("style", "pushbutton", "string", "Update", 'units', 'normalized', "position",[.005 .55-0.3 .12 .05], "callback", @update_cut); 

% Cut column option
h.LCc = uicontrol("style", "text", "string", "Cut cols:", 'units', 'normalized', "position",[.003 .65-0.3 .12 .04]); 
h.PCc = uicontrol("style", "edit", "string", "", 'units', 'normalized', "position",[.005 .60-0.3 .12 .05]); 

%% Colormap selection
%h.LCM = uicontrol("style", "text", "string", "Colormap:", 'units', 'normalized', "position",[.003 .65-0.45 .12 .04]);
%h.PCM = uicontrol("style", "popupmenu", "string", {"rainbow";"hsv"}, 'value', 1, 'units', 'normalized', "position",[.005 .60-0.45 .12 .05]);

% Background color: magenta,(ToDo: cyan, green)
h.LBG = uicontrol("style", "text", "string", "BG color:", 'units', 'normalized', "position",[.003 .06 .08 .04]); 
h.PBG = uicontrol("style", "popupmenu", "string", {'magenta'}, 'value', 1, 'units', 'normalized', "position",[.005 .01 .12 .05]); 

guidata(hf, h);
update_image(hf, true);