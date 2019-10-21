function mamut_xml_from_tracks(tracks, output_path, pixel_size)
% MAMUT_XML_FROM_TRACKS(mamut_file, tracks)
%
%   Convert tracks from a matlab table of tracks to a mamut file.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

if nargin < 3
    pixel_size = [1 1 1];
end

tracks.time = tracks.time - 1;
tracks.x = pixel_size(2)*(tracks.x - 1);
tracks.y = pixel_size(1)*(tracks.y - 1);
tracks.z = pixel_size(3)*(tracks.z - 1);

tracks = tracks{:,:};

[~,d,~] = fileparts(output_path);

% header
doc = com.mathworks.xml.XMLUtils.createDocument('TrackMate');
trackMate = doc.getDocumentElement;
trackMate.setAttribute('version', '2.8.1');
model = doc.createElement('Model');
model.setAttribute('spatialunits', 'pixels');
model.setAttribute('timeunits', 'frames');
trackMate.appendChild(model);

% add spots
allSpots = doc.createElement('AllSpots');
allTimepoints = unique(tracks(:, 8));
for i = 1:numel(allTimepoints)
    spotsInFrame = doc.createElement('SpotsInFrame');
    spotsInFrame.setAttribute('frame', num2str(allTimepoints(i)));
    spots = tracks(tracks(:, 8)==allTimepoints(i), :);
    for j = 1:size(spots, 1)
        spot = doc.createElement('Spot');
        spot.setAttribute('ID', num2str(spots(j, 1)));
        spot.setAttribute('VISIBILITY', '1');
        if spots(j, 6) < 0
            spot.setAttribute('RADIUS', 'NaN');
        else
            spot.setAttribute('RADIUS', num2str(spots(j, 6), '%.1f'));
        end
        spot.setAttribute('QUALITY', num2str(spots(j, 9), '%.1f'));
        spot.setAttribute('POSITION_T', num2str(spots(j, 8), '%.1f'));
        spot.setAttribute('POSITION_X', ...
            num2str(spots(j, 3), '%.13f'));
        spot.setAttribute('POSITION_Y', ...
            num2str(spots(j, 4), '%.13f'));
        spot.setAttribute('FRAME', num2str(spots(j, 8)));
        spot.setAttribute('POSITION_Z', ...
            num2str(spots(j, 5), '%.13f'));
        spotsInFrame.appendChild(spot);
    end
    allSpots.appendChild(spotsInFrame);   
end
model.appendChild(allSpots);

% add tracks
allTracks = doc.createElement('AllTracks');
allTrackID = unique(tracks(:, 10));
for i = 1:numel(allTrackID)
    track = doc.createElement('Track');
    track.setAttribute('name' , ['Track_' num2str(allTrackID(i))]);
    track.setAttribute('TRACK_ID', num2str(allTrackID(i)));
    edges = tracks(tracks(:, 10)==allTrackID(i), :);
    for j = 1:size(edges, 1)
        if edges(j, 7) > 0
            edge = doc.createElement('Edge');
            edge.setAttribute('SPOT_SOURCE_ID', num2str(edges(j, 1)));
            edge.setAttribute('SPOT_TARGET_ID', num2str(edges(j, 7)));
            track.appendChild(edge);
        end
    end
    % if the edge contains a single point, do not add the spot but give a
    % warning instead
    if size(edges, 1) > 1
        allTracks.appendChild(track);
    else
        disp(['Warning: Track ' num2str(allTrackID(i)) ...
            ' contains a single spot, omitted'])
    end
end
model.appendChild(allTracks);

% set tracks
filteredTracks = doc.createElement('FilteredTracks');
for i = 1:numel(allTrackID)
    trackID = doc.createElement('TrackID');
    trackID.setAttribute('TRACK_ID', num2str(allTrackID(i)));
    filteredTracks.appendChild(trackID);
end
model.appendChild(filteredTracks);

% Copy gui state and settings from existing mamut xml file.
mamut_filename = fullfile(output_path, [d '-mamut.xml']);
old_doc = xmlread(mamut_filename);

all_settings = old_doc.getElementsByTagName('Settings');
settings = all_settings.item(0);

all_gui_state = old_doc.getElementsByTagName('GUIState');
gui_state = all_gui_state.item(0);

doc.adoptNode(settings);
doc.adoptNode(gui_state);

trackMate.appendChild(settings.cloneNode(true));
trackMate.appendChild(gui_state.cloneNode(true));

archival_filename = append_timestamp(mamut_filename);
movefile(mamut_filename, archival_filename);
xmlwrite(mamut_filename, doc);
