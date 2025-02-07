function [SCA,T,w_opt,SCM_thresh_opt] = convertRGB2SCA(RGB,datavis,step)
%convertRGB2SCA Uses an input of a Red-Green-Blue (RGB) image and leverages 
%known optical properties of snow to converts the image to a binary snow 
%covered area (SCA) map
%
%   This approach uses information on the optical image color intensity and
%   the differences between RGB bands to detect snow. It relies on the 
%   spectral response of snow which has high reflectance across all optical
%   bands. No color is reflected significantly more than another in most 
%   cases, but shaded snow can present with a blue-ish hue. Using this
%   information, two metrics are computed, combined, then thresholded to
%   produce a snow covered area (SCA, or snow cover) map. The detailed 
%   steps to produce thesese products and outputs are as follows.
%
%
%   Metrics:
%
%   W = normalized distance from white
%   BRd = normalized BLUE - RED band difference
%   w = weighting metric
%   SCM (snow cover metric) = W(w) + BRd(1 - w)
%   SCA (snow covered area or snow cover) = SCM > SCM threshold
%
%
%   Workflow:
%
%   (1) User interactively selects training areas within the larger image
%   and manually tunes the thresholding metrics to produce a reference snow
%   cover map (this is repeated for multiple areas, more = better)
%   
%   (2) These reference maps are used to identify globally optimal values
%   for the metrics 'w' and the 'SCM threshold' by maximize the average
%   F1-score across all delineated reference areas. Optimal values are
%   determined using only the training areas (~75% of delineated areas)
%   while test areas (~25% of delineated areas) are held out for evaluation
%   
%   (3) Error metrics are produced and the optimal thresholds are applied
%   image wide
%
%
%   Inputs:
%
%   RGB (required) - a path or matrix of an image, needs to be 8-bit 
%       (uint8) or 16-bit (uint16) image. Can be georeferenced, or not
%
%   datavis (default: 'on') - if set as 'on', will produce a visualization
%       of the areas delineated as snow covered by the function
%
%   step (default: 0.05) - sets the step size to use in the optimization
%       procedure. Must be in the range 0.01 (very slow) to 0.1. The 
%       default value works well in most cases. If results are poor using
%       the default, try decreasing towards 0.01 for improvement
%
%
%   Outputs:
%
%   SCA - the produced 2-D snow cover grid as numeric uint8, where: 
%       0 = no snow, 1 = snow, 255 = unclassified/fill areas. This product 
%       is also saved at the same path as the input image with the 
%       additional suffix '_SCA' as a .tif image. If the input is 
%       georeferenced, this output will be as well
%
%   T - a table containing performance metrics and characteristics of each
%       reference area (test and training) as well as the averages. This is 
%       saved to the same path as the input image with the additional 
%       suffix '_SCA_evaluation' as a .csv file
%
%   w_opt - the optimal weighting metric (w) value used to produce 'SCA'
%
%   SCM_thresh_opt - the optimal SCM threshold value used to produce 'SCA'
%
%
% Notes:
%
% 1) tested on images created from a Phantom 4 Pro and Iphone 11 camera. 
%    Band center wavelengths are ~594nm (570-680nm, Red), ~525nm 
%    (470-620nm, Green), and ~460nm (400-550nm, Blue) [from Burggraaff et 
%    al., 2019 for P4 Pro Camera]
%
% 2) a diverse range of test areas leads to a better SCA result for the
%    full image
%
% 3) quality of the output SCA image is dependent on the analysts ability
%    to effectively delineate snow covered areas over test areas
%
% 4) the interactive thresholding functionality may not work for all
%    environments, as it has been minimally tested in scenes including both
%    clouds and snow cover and was designed to classify ground facing 
%    (e.g., UAS survey) imagery. The method is challenged if snow, shaded 
%    snow, clouds, and blue sky are present
%
% 5) the function requires two additional non-built in dependencies: 
%    (1) sliderClassifySCA and (2) F1score. These must also be added to the
%    user's path to use the convertRGB2SCA function
%
% 6) program is memory intensive, and may require sub-setting for large 
%    images (generally >20,000 x 20,000). However, this is computer/MATLAB 
%    settings dependent

%{
%use for testing
RGB = '/Users/jjohns/Pictures/DSCF1608.JPG';
step = 0.05; 
datavis = "on";
%}

%set defaults
narginchk(1,3);
if nargin == 1
    datavis = "on";
    step = 0.05;
elseif nargin == 2
    step = 0.05;
end

%constrain inputs to valid range
if step < 0.01
    step = 0.01;
elseif step > 0.1
    step = 1;
end

%% (1) Load in imagery and identify areas to mask
tic;
disp(' ')
disp('Loading imagery...')

%determine if is path to a georeferenced image, non-georeferenced image. 
% If not, is considered a numeric matrix input
if ischar(RGB) | isstring(RGB)
    %input is a path to file
    RGB_file = RGB;

    %attempt to load geotiff
    try
        [RGB,REF] = readgeoraster(RGB_file);

        %get the EPSG code from the user, as it is not included in the
        %MapCellsReference object
        try
            EPSGcodeName = REF.ProjectedCRS.Name;
            prompt = {strcat("Enter EPSG code for ",EPSGcodeName)};
            dlgTitle = 'User Input';
            fieldsize = [1 45];
            EPSG = inputdlg(prompt,dlgTitle,fieldsize);
            EPSG = str2double(EPSG{:});
            isGeographic = 0;
        if isnan(EPSG)
            error('Invalid EPSG code')
        end
        %if is in geographic coordinates, will save with REF object
        catch
            isGeographic = 1;
        end
        
        
    catch
        %read in as normal image
        RGB = imread(RGB_file);
    end
else
    %check that is at least a 3-band image
    if size(RGB,3) < 3
        error('Input is not a 3-band image');
    end
end

%determine RGB colorspace (only 8-bit or 16-bit rgb right now)
if isa(RGB,"uint8")
    M = 2^8 - 1;
elseif isa(RGB,"uint16")
    M = 2^16 - 1;
else
    error('Unsupported filetype')
end

%determine masked areas (dark or as defined in image alpha band)
if size(RGB,3) == 4
    mask = RGB(:,:,4) == 0;
else 
    mask = zeros(size(RGB(:,:,1),[1 2]),'logical') ~= 0;
end    

%display progress
n = toc;
disp(['Image load completed in ' num2str(n) ' seconds']);


%% (2) Define SCA reference areas
tic;
disp(' ')
disp('Defining reference SCA areas...')

%visualize the image
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(RGB(:,:,1:3),'AlphaData',~mask);
title('Click within the image and drag to delineate a reference area','Fontsize',18)

%loop through creating test areas until user terminates the loop
test_areas = 0;
while true
    %keep track of iterations through the loop
    test_areas = test_areas + 1;

    %user selects an ROI
    roi = drawrectangle;
    pos = round(roi.Position); %get rounded position information

    %define bounds in terms of row/col positions
    row_min = pos(2);
    row_max = pos(2)+pos(4);
    col_min = pos(1);
    col_max = pos(1)+pos(3);

    %get cropped version of the image based on user input
    RGBc = RGB(row_min:row_max,col_min:col_max,1:3);

    %do interactive SCA mapping and return local metrics
    h = sliderClassifySCA(RGBc);

    %add test area positions to table & store the cropped images in a
    %dynamically updated structure
    if test_areas == 1        
        %create and update structure with test area data
        im_struct = struct;
        im_struct.RGB.(['ID' num2str(test_areas)]) = h.RGBc;
        im_struct.W.(['ID' num2str(test_areas)]) = h.W;
        im_struct.BRd.(['ID' num2str(test_areas)]) = h.BRd;
        im_struct.SCA.(['ID' num2str(test_areas)]) = h.SCA;
        im_struct.bounds.(['ID' num2str(test_areas)]) = [row_min row_max col_min col_max];       
    else  
        %update structure with test area data
        im_struct.RGB.(['ID' num2str(test_areas)]) = h.RGBc;
        im_struct.W.(['ID' num2str(test_areas)]) = h.W;
        im_struct.BRd.(['ID' num2str(test_areas)]) = h.BRd;
        im_struct.SCA.(['ID' num2str(test_areas)]) = h.SCA;
        im_struct.bounds.(['ID' num2str(test_areas)]) = [row_min row_max col_min col_max];
    end

    %repeat until user does not want to define another ROI
    dlgTitle    = 'User Question';
    dlgQuestion = 'Select another test area?';
    choice = questdlg(dlgQuestion,dlgTitle,'Yes','No, classify image', 'Yes');
    if strcmp(choice,'No, classify image')
        break
    end
end

%remove figure showing RGB image & test areas
delete(gcf);

%display progress
n = toc;
disp(['Reference areas created in ' num2str(n) ' seconds']);


%% (3) Identify the optimal weighting value (w) and SCM threshold
tic;
disp(' ')
disp('Running threshold optimization...')

%check that enough areas are selected
if test_areas < 2
    error('Not enough training areas: user must delineate at least 2 training areas')
end

%partition the reference areas into training and test areas
N_holdout = round(test_areas .* 0.25);
if N_holdout < 1
    N_holdout = 1;
end
IDs_testing = sort(randperm(test_areas,N_holdout));
IDs_training = 1:test_areas;
IDs_training = setxor(IDs_training,IDs_testing);

%store iteration counts
iter = 0;
%loop through valid range of weighting metrics at intervals of size step
for w = 0:step:1
    %loop though a range of SCM thresholds from -0.1 (BRd can be negative 
    % in rare cases when the red channel intensity > blue) to 1 (max)
    for SCM_thresh = -0.1:step:1

        %iteration counter
        iter = iter + 1;

        %compute and store classification quality metrics
        F1_scores = NaN([test_areas 1]);
        precision_scores = NaN([test_areas 1]);
        recall_scores = NaN([test_areas 1]);
        accuracy_scores = NaN([test_areas 1]);

        %iterate through each training area
        for i = 1:length(IDs_training)
            
            %get id of training area
            ID_training = IDs_training(i);

            %get ID for extracting variables from structure
            ID = ['ID' num2str(ID_training)];

            %get reference SCA
            SCA_ref = im_struct.SCA.(ID);

            %calculate SCA using thresholds for test area i
            SCA_i = ((im_struct.W.(ID) .* w) + (im_struct.BRd.(ID) .* (1 - w))) > SCM_thresh;

            %calculate the F1 score comparing the reference classification
            %to the thresholded one (SCA_i)
            [F1,precision,recall,accuracy] = F1score(SCA_ref,SCA_i);

            %update arrays with the computed scores
            F1_scores(i) = F1;
            precision_scores(i) = precision;
            recall_scores(i) = recall;
            accuracy_scores(i) = accuracy;

        end

        %update results table
        if iter == 1

            %create table
            T = table();
            T.w = w;
            T.SCM_thresh = SCM_thresh;          
            T.F1_mean = mean(F1_scores,"all","omitmissing");

        else
            
            %update table with new row
            Ti = table();
            Ti.w = w;
            Ti.SCM_thresh = SCM_thresh;
            Ti.F1_mean = mean(F1_scores,"all","omitmissing");
            T = [T; Ti];

        end          
    end
end

%Generate a 2D interpolated field at higher resolution using the spline
%technique. Spline is better for extrapolating between points and results
%in improved estimates of the global maxima. For example, linear
%interpolation does not result in unique values between points and the
%optimal solution must therefor be one from the prior step, making the 
%interpolation step arbitrary. This step saves substantial computational
%time that would be required if calculating results of all metric
%combinations directly

%generate a target grid with high spatial density of points (step of 0.001)
w_vals = 0:0.001:1;
SCM_thresh_vals = -0.1:0.001:1;
[xq,yq] = meshgrid(w_vals,SCM_thresh_vals);

%input the tested weighting values, thresholds, and F1 scores
x = T.w;
y = T.SCM_thresh;
v = T.F1_mean;
idx = ~isnan(v);
vq = griddata(x(idx),y(idx),v(idx),xq,yq,'v4');

%identify the location of the global maxima
[~,ind] = max(vq(:));
[r,c] = ind2sub(size(vq),ind);

%get metric values at the global maxima (optimal metric values)
w_opt = xq(r,c);
SCM_thresh_opt = yq(r,c);

%visualize the result
figure;
%plot all tested points
scatter3(x,y,v,[],v,'filled','HandleVisibility','off');
xlabel('weighting value (w)');
ylabel('SCM threshold','Interpreter','none');
zlabel('mean F1-score');
cb = colorbar;
ylabel(cb,'mean F1-score');
hold on; 
%add global maxima to plot
scatter3(w_opt,SCM_thresh_opt,vq(r,c),75,'r','filled');
set(gca,'Fontsize',14);
title(['Point of Optimality (red): w = ' num2str(w_opt,3) ', SCM threshold =' num2str(SCM_thresh_opt,3)]);

%display progress
n = toc;
disp(['Threshold optimization completed in ' num2str(n) ' seconds']);


%% (4) Apply optimized metric values to create and save final SCA products
tic
disp(' ')
disp('Classifying image and saving snow cover products...')

%calculate metrics from RGB image
%get individual normalized RGB bands for valid areas
R = single(RGB(:,:,1))./M; R = R(~mask);
G = single(RGB(:,:,2))./M; G = G(~mask);
B = single(RGB(:,:,3))./M; B = B(~mask);

%calculate normalized euclidean distance from ideal white (W)
W = (((1 - R).^2) + ((1 - G).^2) + ((1 - B).^2)).^0.5;
W = (W./((3*(1.^2)).^0.5)); %normalize (0-1) based on max possible distance
W = (1 - W); %invert so higher W indicates whiter surfaces in the image

%calculate blueness metric as the normalized difference between the 
%blue and red bands as shaded snow presents with a blue-ish hue
BRd = (B - R);

%clear variables to free up memory
clear R G B;

%convert metrics to grids
W_grid = NaN(size(mask),'single');
W_grid(~mask) = W;

%convert metric to grid
BRd_grid = NaN(size(mask),'single');
BRd_grid(~mask) = BRd;

%clear metrics to free up memory
clear W BRd;

%apply optimized thresholds to the full image to create SCA map
SCA = ((W_grid .* w_opt) + (BRd_grid .* (1 - w_opt))) > SCM_thresh_opt;
SCA = uint8(SCA); %convert to uint8
SCA(mask) = 255; %assign 255 to invalid/unclassified areas

%save the SCA map as a georeferenced or un-georeferenced image, depending
%on input type
if exist("REF","var")
    if isGeographic == 1
        geotiffwrite([RGB_file(1:end-4) '_SCA.tif'],SCA,REF);
    elseif isGeographic == 0
        geotiffwrite([RGB_file(1:end-4) '_SCA.tif'],SCA,REF,"CoordRefSysCode",EPSG);
    end
else
    imwrite(SCA,[RGB_file(1:end-4) '_SCA.tif']);
end

%compute and store classification quality metrics for optimal case
bounds = cell([test_areas 1]);
px_counts = NaN([test_areas 1]);
test_area_ids = NaN([test_areas 1]);
F1_scores = NaN([test_areas 1]);
precision_scores = NaN([test_areas 1]);
recall_scores = NaN([test_areas 1]);
accuracy_scores = NaN([test_areas 1]);

%loop through test areas and recalculate performance
for i = 1:test_areas

    %store information on if the area was used for training or testing,
    %test areas will be given a 1, training areas used in optimization are
    %given a 0
    if sum(i == IDs_testing) == 1
        test_area_ids(i) = 1;
    else
        test_area_ids(i) = 0;
    end

    %get ID for extracting variables from structure
    ID = ['ID' num2str(i)];

    %get reference SCA
    SCA_ref = im_struct.SCA.(ID);

    %calculate SCA using thresholds for test area i
    SCA_i = ((im_struct.W.(ID) .* w_opt) + (im_struct.BRd.(ID) .* (1 - w_opt))) > SCM_thresh_opt;

    %calculate the F1 score comparing the reference classification
    %to the thresholded one (SCA_i)
    [F1,precision,recall,accuracy] = F1score(SCA_ref,SCA_i);

    %update arrays with the computed scores
    F1_scores(i) = F1;
    precision_scores(i) = precision;
    recall_scores(i) = recall;
    accuracy_scores(i) = accuracy;
    px_counts(i) = numel(SCA_ref);
    bounds{i} = im_struct.bounds.(ID);

end

%create table with all performance information
T = table();
T.w_opt = repmat(w_opt,size(F1_scores));
T.SCM_thresh_opt = repmat(SCM_thresh_opt,size(F1_scores));
T.test_area = test_area_ids;
%assign each cell bounding position a table column [row_min row_max col_min col_max]
fun = @(x) x(1);
T.row_min = cellfun(fun,bounds);
fun = @(x) x(2);
T.row_max = cellfun(fun,bounds);
fun = @(x) x(3);
T.col_min = cellfun(fun,bounds);
fun = @(x) x(4);
T.col_max = cellfun(fun,bounds);
T.pixel_counts = px_counts;
T.precision_scores = precision_scores;
T.recall_scores = recall_scores;
T.accuracy_scores = accuracy_scores;
T.F1_scores = F1_scores;

%append a row containing averages
Ti = table();
Ti.w_opt = NaN;
Ti.SCM_thresh_opt = NaN;
Ti.test_area = NaN;
Ti.row_min = NaN;
Ti.row_max = NaN;
Ti.col_min = NaN;
Ti.col_max = NaN;
Ti.pixel_counts = mean(px_counts,"all","omitmissing");
Ti.precision_scores = mean(precision_scores,"all","omitmissing");
Ti.recall_scores = mean(recall_scores,"all","omitmissing");
Ti.accuracy_scores = mean(accuracy_scores,"all","omitmissing");
Ti.F1_scores = mean(F1_scores,"all","omitmissing");
T = [T; Ti];

%return SCA product quality evaluation .csv file
writetable(T,[RGB_file(1:end-4) '_SCA_evaluation.csv'],'WriteMode','overwrite');

n = toc;
disp(['Products created and saved in ' num2str(n) ' seconds']);


%% (5) Optional: Visualize final SCA product
%create visualization showing snowy areas
if strcmp(datavis,"on")
    pause(1)
    %visualize the RGB image with snow covered areas removed
    figure;
    %normal image, but resized to having only 2000 rows (for display)
    Nrow = 2000;
    Ncol = round(size(SCA,2)/(size(SCA,1)/Nrow));
    imshow(imresize(RGB(:,:,1:3),[Nrow Ncol]));
    %overlay mask of snow covered areas (in blue) and invalid areas (red)
    alphamask(imresize(SCA == 1,[Nrow Ncol]), [0 0 1], 0.5); hold on;
    alphamask(imresize(SCA == 255,[Nrow Ncol]), [1 0 0], 0.5);
    title(['Blue = snow cover, Red = masked areas   ',...
        '(w: ' num2str(w_opt,3) ', SCM threshold: ' num2str(SCM_thresh_opt,3)...
        ', F1-score: ' num2str(T.F1_scores(end),3) ')'],...
        'Interpreter','none','FontSize',14);
end

end