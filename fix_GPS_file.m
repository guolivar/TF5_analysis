clear all
fclose all
# close all
page_screen_output(0);
page_output_immediately(1);
# Set the scene
work_dir='/home/olivaresga/data/TF5_JapanNZvoyage/2013_May/';
work_file='GPS_201305.txt';
fixed_file=[work_file(1:end-4) '_dat.txt'];
header_file=[work_file(1:end-4) '_hdr.txt'];
id_in=fopen([work_dir work_file],'rt');
id_ou=fopen([work_dir fixed_file],'w');
id_hdr=fopen([work_dir header_file],'w');
# Read lines, send the NON DATA to header as they are and parse the DATA lines
# DATA lines: line(1)=='F' | line(1)=='T'
# NON DATA: everything else
while ~feof(id_in)
  line=fgetl(id_in);
  if length(line)>0
    if (line(1)=='F' | line(1)=='T')
      if line(1)=='T'
	idx=strfind(line,'S');
	if (length(idx)>0)
	  line(idx(1))='-';
	end
	idx=strfind(line,'N');
	if (length(idx)>0)
	  line(idx(1))=' ';
	end
	idx=strfind(line,'W');
	if (length(idx)>0)
	  line(idx(1))='-';
	end
	idx=strfind(line,'E');
	if (length(idx)>0)
	  line(idx(1))=' ';
	end
	fprintf(id_ou,'%s\n',line);
      else
	fprintf(id_ou,'%s\n',line);
      end
    else
      fprintf(id_hdr,'%s\n',line);
    end
  end
end
fclose all