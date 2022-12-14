%% Cell segmentation and reads assignment based on overlap
% prototype
% Cartana, Xiaoyan, 2021-4-12

nuclei_image = 'R:\Benutzer\Sallinger Katja\Shedder project\shedder\260KS\260KS_AP1_DAPI.tif';
reads_file = 'R:\Benutzer\Sallinger Katja\Shedder project\shedder\260KS\Decoded_LowThreshold.csv';
output_prefix = 'R:\Benutzer\Sallinger Katja\Shedder project\shedder\260KS\Segmentation\output';

%% load image and coordiantes
DAPI = imread(nuclei_image);
spots = importdata(reads_file);
gene = spots.textdata(2:end,1);
xy = spots.data;

if size(DAPI,3) > 1
    DAPI = rgb2gray(DAPI);
end

% scenario 1 - enough RAM to process the full image
% 
% binary thresholding (it works generally well without smoothing)
% Otsu's method, threshold decided by the full DAPI image
% Ibw = im2bw(DAPI, graythresh(DAPI));
% 
% h-minima filter to suppress pixels that have less than 0.001
% (converted to 0-1 range) depth, after smoothing using a Guassian
% filter (generally to avoid oversegmentation in watershed)
% if isa(DAPI, 'uint16')
%     Ismooth = imfilter(double(DAPI)/65535, fspecial('gaussian', 10, 2));
%     Ihmin = imhmin(-Ismooth, .001);
% else
%     Ismooth = imfilter(double(DAPI)/255, fspecial('gaussian', 10, 2));
%     Ihmin = imhmin(-Ismooth, .1);
% end    
% 
% watershed on h-minima filetered image
% Iws = watershed(Ihmin);
% 
% impose the binay mask created from thresholding
% Iws = double(Iws) .* double(Ibw);
% 
% expand 30 px to mimic cells
% [D, idx] = bwdist(Iws);
% D = uint32(D<=30);
% Icell = reshape(Iws(idx(:)), size(Iws));
% Icell = D.*uint32(Icell);
% 
% get centroid and area of each cell
% props = regionprops(Icell, 'Centroid', 'Area');
% Cells =  [(1:length(props))', cat(1, props.Centroid), cat(1, props.Area)];
% 
% outline image
% Outline = Icell ~= imerode(Icell, strel('disk',1));
% 
% determine spot-cell relation based on overlap
% Parent = Icell(sub2ind([size(DAPI,1), size(DAPI,2)], round(xy(:,2)), round(xy(:,1))));

%% scenario 2 - not enough RAM to process the whole image
% major steps are the same as above

% determine how many tiles in X and Y direction
nX = ceil(size(DAPI,2)/2000);
nY = ceil(size(DAPI,1)/2000);

TileSpotCell = cell(nX*nY,1);
% EdgeObjects = cell(nX*nY,4);
Boundaries = cell(nX*nY,1);
CellProps = cell(nX*nY,3);

parfor t = 1: (nX*nY)
    fprintf('Processing tile %d/%d\n', t, nX*nY);
    
    [iX, iY] = ind2sub([nX, nY], t);    % counting of t starts with x axis
    
    % tile images into 2000 +- 200 pixel size
    % extra 200 is for easy edge correction
    I = DAPI(max(1, 2000*(iY-1)-200+1) : min(size(DAPI,1), 2000*iY+200),...
        max(1, 2000*(iX-1)-200+1) : min(size(DAPI,2), 2000*iX+200));
    padX0 = 200 * (1 < 2000*(iX-1)-200+1);
    padY0 = 200 * (1 < 2000*(iY-1)-200+1);
    padX1 = 200 * (size(DAPI,2) > 2000*iX+200);
    padY1 = 200 * (size(DAPI,1) > 2000*iY+200);
    
    szI = size(I);
    padX1 = max(szI(2)-2200, padX1);
    padY1 = max(szI(1)-2200, padY1);
        
    % binary thresholding (it works generally well without smoothing)
    % Otsu's method, threshold decided by the full DAPI image
    Ibw = im2bw(I, graythresh(DAPI));

    if nnz(Ibw)
        % h-minima filter to suppress pixels that have less than 0.0001 
        % (16-bit converted to 0-1 range) depth after smoothing using a 
        % Guassian filter (generally to avoid oversegmentation in watershed)
        Ismooth = imfilter(double(I)/65535, fspecial('gaussian', 10, 2));
        Ihmin = imhmin(-Ismooth, .0001);

        % watershed on h-minima filetered image
        Iws = watershed(Ihmin);

        % impose the binay mask created from thresholding, ignore all
        % background pixels
        Iws = double(Iws) .* double(Ibw);

        % expand 30 px to mimic cells
        [D, idx] = bwdist(Iws);
        D = uint32(D<=30);
        Icell = reshape(Iws(idx(:)), size(Iws));
        Icell = D.*uint32(Icell);
        
        % get centroid and area of each cell
        props = regionprops(Icell, 'Centroid', 'Area');
        
        % remove any cell whose centroid is not in the center 2000x2000 px
        inCenter = all(cat(1,props.Centroid) > [padX0, padY0],2)...
            & all(cat(1,props.Centroid) <= fliplr(szI) -[padX1, padY1],2);
        
        % only cells in the center are left
        if nnz(inCenter)        
            Icell(ismember(Icell, find(~inCenter))) = 0;
            
            CellProps(t,:) = [{cat(1, props(inCenter).Centroid)}, {cat(1, props(inCenter).Area)}, {find(inCenter)}];

            % assign reads to cells based on overlap    
            xyTile = round(xy) - [max(0, 2000*(iX-1)-200), max(0, 2000*(iY-1)-200)];
            inTile = all(xyTile>0,2) & xyTile(:,1)<=szI(2) & xyTile(:,2)<=szI(1);
            parent = Icell(sub2ind(szI, xyTile(inTile,2), xyTile(inTile,1)));

            TileSpotCell{t} = [find(inTile), parent];

    %         % DAPI cut by north, west, south and east lines
    %         EdgeObjects(t,:) = [{[find(Iws(1,:))', Iws(1,Iws(1,:)~=0)']},...
    %             {[find(Iws(:,1)), Iws(Iws(:,1)~=0,1)]},...
    %             {[find(Iws(szI(1),:))', Iws(szI(1),Iws(szI(1),:)~=0)']},...
    %             {[find(Iws(:,szI(2))), Iws(Iws(:,szI(2))~=0,szI(2))]}];

            % save boundary image
            Ibound = Icell ~= imerode(Icell, strel('disk',1));
            Boundaries{t} = Ibound(padY0+1 : szI(1)-padY1, padX0+1 : szI(2)-padX1);
        end
    end
end

% assemble the outline image
Outline = false(size(DAPI,1), size(DAPI,2));
for t = 1:(nX*nY)
    [iX, iY] = ind2sub([nX, nY], t);    % counting of t starts with x axis
    
    if ~isempty(Boundaries{t})
        Outline(max(1, 2000*(iY-1)+1) : min(size(DAPI,1), 2000*iY),...
        max(1, 2000*(iX-1)+1) : min(size(DAPI,2), 2000*iX)) = Boundaries{t};
    end
end

% renumber and organize output info
Cells = [];
IdxConversion = [];
nCell = 0;
for t = 1:(nX*nY)
    if ~isempty(CellProps{t,1})
        [iX, iY] = ind2sub([nX, nY], t);    % counting of t starts with x axis
        
        cellpos = CellProps{t,1};
        cellarea = CellProps{t,2};
        
        % conversion between TileSpotCell and final Cell index
        IdxConversion = [IdxConversion;...
            repmat(t, size(cellpos,1),1), CellProps{t,3}, nCell+(1:length(CellProps{t,3}))'];
        
        % global position of cells
        Cells = [Cells;...
            nCell+(1:length(CellProps{t,3}))',...
            cellpos + [max(1, 2000*(iX-1)-200+1),max(1, 2000*(iY-1)-200+1)],...
            cellarea];
        nCell = nCell + length(CellProps{t,3});
    end
end

% finalize spot-cell assignment
% if a spot is assigned to more than one cell, the later assignment has
% absolute priority unless it's assigned to 'Cell 0'
Parent = zeros(length(xy),1);
for t = 1:(nX*nY)
    if ~isempty(TileSpotCell{t})
        tempConversion = IdxConversion(IdxConversion(:,1)==t,2:3);
        idx = cellfun(@(v) find(v==tempConversion(:,1)), num2cell(TileSpotCell{t}(TileSpotCell{t}(:,2)~=0,2)));
        Parent(TileSpotCell{t}(TileSpotCell{t}(:,2)~=0,1)) = tempConversion(idx,2);
    end
end

%% in any case, output
% cell properties
fid = fopen([output_prefix '_Cells.csv'], 'w');
fprintf(fid, 'cellID,X,Y,area\n');
fprintf(fid, '%d,%f,%f,%f\n', Cells');
fclose(fid);
    
% spot with parent cell
fid = fopen([output_prefix '_SpotWithCell.csv'], 'w');
fprintf(fid, 'gene,X,Y,cellID\n');
towrite = [gene, num2cell(xy), num2cell(Parent)]';
fprintf(fid, '%s,%f,%f,%d\n', towrite{:});
fclose(fid);

% gene x cell matrix
[uGenes, ~, iGene] = unique(gene);
[uParents, ~, iParent] = unique(Parent);
count = hist3([iParent, iGene], [{1:length(uParents)}, {1:length(uGenes)}]);
fid = fopen([output_prefix '_GeneXCell.csv'], 'w');
header = [{'cellID'}; uGenes; {'cellX' 'cellY' 'cellArea'}']; 
fprintf(fid, ['%s' repmat(',%s', 1, numel(header)-1) '\n'], header{:});
if uParents(1) == 0
    towrite = [uParents, count, [[0 0 0]; Cells(uParents(2:end),2:4)]]';
else
    towrite = [uParents, count, Cells(uParents,2:4)]';
end
fprintf(fid, ['%d' repmat(',%d', 1, numel(uGenes)) ',%f,%f,%d\n'], towrite);
fclose(fid);

% outline image
imwrite(Outline, [output_prefix '_Outline.jpg']);

% outline on dapi
if isa(DAPI, 'uint16')
    imwrite(uint8(DAPI/250) + repmat(uint8(Outline)*255,1,1,3), [output_prefix '_OutlineOnDAPI.jpg']);
else
    imwrite(DAPI + repmat(uint8(Outline)*255,1,1,3), [output_prefix '_OutlineOnDAPI.jpg']);
end
    
    