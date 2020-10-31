function [data] = read_STK(filename)
%   Read and parse an STK ECEF satellite position file (assumes data directory)
%   ----------------------------------------------------------------------------
%   
%   Inputs:
%   --------
%   'filename' -- string of the file name
%     
%   Returns:
%   --------
%   'data' -- struct containing the STK file info

% month container (not used)
% month = containers.Map({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',...
%     'Sep','Oct','Nov','Dec'},{1,2,3,4,5,6,7,8,9,10,11,12});

fid = fopen(filename);

line = fgetl(fid);
i = 0;
j = 0;
while ischar(line)
    if (i >= 7)
        j = j + 1;
        Line = split(line);
        
        % time stuff
        t = [char(Line(1)), char(Line(2)), char(Line(3)), char(Line(4))];
        
        % space stuff
        x = str2double(Line(5));
        y = str2double(Line(6));
        z = str2double(Line(7));
        r_vec = [x, y, z];
        r_mag = str2double(Line(8));

        % append time and space stuff to data
        data.Time(j) = datetime(t, 'InputFormat', 'ddMMMyyyyHH:mm:ss.SSS');
        data.r_vec(j,:) = r_vec;
        data.r_mag(j) = r_mag;
    end
    i = i + 1;
    line = fgetl(fid);
end
fclose(fid);

end

