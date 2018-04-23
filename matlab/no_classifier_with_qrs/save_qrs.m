function save_qrs(data_full,fs,method)
annot = qrs_detector_wrap(data_full,fs,method);
if isequal(method,'gqrs') || isequal(method,'sqrs') || isequal(method,'wqrs')
    annot = qrs_corrector_wrap(data_full,annot,fs);
end

% save output
[pathstr,filename,~] = fileparts(data_full);
outfilename = [pathstr '\' filename '.adj'];
fp = fopen(outfilename, 'wt');
for rPeakIdx=1:length(annot)
    fprintf(fp,'%ld \n', round(annot(rPeakIdx)));
end
fclose(fp);