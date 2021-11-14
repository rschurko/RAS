% function fid = bread(fname)
%	bread.m - read data from Bruker fid or ser files
%		fname:	'fid' or 'ser'
%		fid:	data returned 

function fid = bread(fname);

f=fopen(fname,'r','ieee-le'); 
fid=fread(f,'int32');
td=size(fid,1);
fid=reshape(fid,2,td/2);
fid=fid(1,:)+i*fid(2,:);
fid=reshape(fid,td/2,1);
fclose(f);

