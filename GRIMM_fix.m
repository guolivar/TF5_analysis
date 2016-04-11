clear all
# fclose all
# close all
# Read the pasted GRIMM file
page_screen_output(0);
page_output_immediately(1);
# Set the scene
work_dir='/home/olivaresga/data/TF5_JapanNZvoyage/2013_May/';
work_file='Grimm_201305.TXT';
fixed_file=[work_file(1:end-4) '_fix.txt'];
id_in=fopen([work_dir work_file],'rt');
id_ou=fopen([work_dir fixed_file],'w');
# Write headers
first_line={'OctaveDate','Year','Month','Day','Hour','Minute','Error', ...
  '265', ...
'290', ...
'324', ...
'374', ...
'424', ...
'474', ...
'539', ...
'614', ...
'675', ...
'748', ...
'894', ...
'1140', ...
'1442', ...
'1789', ...
'2236', ...
'2739', ...
'3240', ...
'3742', ...
'4472', ...
'5701', ...
'6982', ...
'7984', ...
'9220', ...
'11180', ...
'13693', ...
'16202', ...
'18708', ...
'22361', ...
'27386', ...
'30984', ...
'35000'};
dp=[265 290 324 374 424 474 539 614 675 748 894 1140  1442  1789  2236  2739  3240  3742  4472  5701  6982  7984  9220  11180 13693 16202 18708 22361 27386 30984 35000];
dlogDp=[0.0492180227  0.0299632234  0.0669467896  0.057991947 0.0511525224  0.0457574906  0.0644579892  0.0494853631  0.0321846834  0.057991947 0.096910013 0.1139433523  0.0901766303  0.096910013 0.096910013 0.079181246 0.0669467896  0.057991947 0.096910013 0.1139433523  0.0621479067  0.0543576623  0.0705810743  0.096910013 0.079181246 0.0669467896  0.057991947 0.096910013 0.079181246 0.0280287236  1];
fmt_hdr='%s';
for i=1:36
  fmt_hdr=[fmt_hdr '\t%s'];
end
fmt_hdr=[fmt_hdr '\t%s\n'];
fmt_data='%f\t%d\t%d\t%d\t%d\t%d\t%d';
for i=1:30
  fmt_data=[fmt_data '\t%f'];
end
fmt_data=[fmt_data '\t%f\n'];
pad_data=NaN*ones(1,31);
fprintf(id_ou,fmt_hdr,first_line{:});
# Read the first lines until encountering a P line
c_line=fgetl(id_in);
while and(~feof(id_in),c_line(1)!='P')
  c_line=fgetl(id_in);
  c_line
end
while ~feof(id_in)
  #The starting point is a P line that's already in c_line
  #Get date, time and error code
  p_vec=sscanf(c_line(3:end),'%f')';
  if length(p_vec)<10
    p_vec=NaN*ones(1,10);
    octave_date=NaN;
  else
    octave_date=datenum([p_vec(1:5) 0]);
  end
  data_vec1=[octave_date p_vec(1:5) p_vec(8)];
  #Get the number concentration lines (C records)for the whole minute
  next_P=(0==1);
  invalid=0;
  data_vec2=zeros(1,31);
  for imin=1:1
    #get 4 lines of data and populate the data vector
    #1
    if ~feof(id_in)
      c_line=fgetl(id_in);
    else
      invalid=1;
      break;
    end
    if and(c_line(1)=='C',length(c_line)>73)
      data_vec2(1:8)=data_vec2(1:8)+0.1*sscanf(c_line(4:74),'%f')';
    else
      invalid=1;
      break;
    end
    #2
    if ~feof(id_in)
      c_line=fgetl(id_in);
    else
      invalid=1;
      break;
    end
    if and(c_line(1)=='C',length(c_line)>73)
      #This line has the repeated 2500nm channel
      data_vec2(9:15)=data_vec2(9:15)+0.1*sscanf(c_line(4:65),'%f')';
    else
      invalid=1;
      break;
    end
    #3
    if ~feof(id_in)
      c_line=fgetl(id_in);
    else
      invalid=1;
      break;
    end
    if and(c_line(1)=='c',length(c_line)>73)
      data_vec2(16:23)=data_vec2(16:23)+0.1*sscanf(c_line(4:74),'%f')';
    else
      invalid=1;
      break;
    end
    #4
    if ~feof(id_in)
      c_line=fgetl(id_in);
    else
      invalid=1;
      break;
    end
    if and(c_line(1)=='c',length(c_line)>73)
      data_vec2(24:31)=data_vec2(24:31)+0.1*sscanf(c_line(4:74),'%f')';
    else
      invalid=1;
      break;
    end
  end
  if ~isnan(data_vec1(1))
    #The DATE is valid (full P record) process
    if invalid
      #There were not 10x4 lines of data, i.e., not a full minute
      #So, print an invalid record
      fprintf(id_ou,fmt_data,[data_vec1 pad_data]);
    else
      #We have a full minute so print the record to file
      #First, we calculate the dN/dlogDp for each channel,
      #except the last one that keeps as is
      data_vec2(1:30)=data_vec2(1:30)-data_vec2(2:31);
      data_vec2=data_vec2./dlogDp;
      fprintf(id_ou,fmt_data,[data_vec1 data_vec2]);
    end
  else
    #The DATE is invalid (not full P record) don't do anything
  end
  #Now look for the next P line
  next_P=0;
  while and(~feof(id_in),~next_P)
    c_line=fgetl(id_in);
    if length(c_line)>75
      next_P=c_line(1)=='P';
    else
      next_P==0;
    end
  end
end
fclose(id_in);
fclose(id_ou);
