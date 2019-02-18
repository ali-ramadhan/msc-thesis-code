function M = geometryFilesToMatrix(filenamePrefix, indices, debug)
  nFiles = size(indices,2);
  M = [];

  for i = indices
    if (debug && (rem(i, 100) == 0)) fprintf('Reading file %d/%d.\n', i, size(indices,2)); end;

    filename = strcat(filenamePrefix, '_G', sprintf('%05d',i), '.log');
    fIn = fopen(filename, 'r');

    formatSpec = '%d %f %f %f %f %d';
    fileSize = [6 Inf];
    geometries = fscanf(fIn, formatSpec, fileSize);

    fclose(fIn);

    M = [M; (geometries)';]; % TODO: This is inefficient but probably okay for us.
  end

  nGeometries = size(M,1);
  fprintf('Imported %d geometries from %d individual geometry files.\n', nGeometries, nFiles);
  fprintf('Don''t forget to filter out bad geometries (e.g. exitFlag != 2).\n');

  prompt = 'Save all geometries into one file? Y/N [Y]: ';
  saveOneFile = input(prompt, 's');
  if isempty(saveOneFile) str = 'Y'; end;

  if strcmp(saveOneFile, 'Y')
    defaultFilename = strcat(filenamePrefix, '.txt');
    prompt = sprintf('Enter file name [%s]: ', defaultFilename);
    fOutFilename = input(prompt, 's');
    if isempty(fOutFilename) fOutFilename = defaultFilename; end;

    fprintf('Saving %d individual geometry files to %s (%d geometries contained)...\n', nFiles, fOutFilename, nGeometries)

    fOut = fopen(fOutFilename, 'a');
    for i = 1:nGeometries
      fprintf(fOut, '%d\t%3.6f\t%3.6f\t%3.6f\t%2.2f\t%d\n', M(i,:));
    end
    fclose(fOut);
  end

  prompt = 'Delete individual geometry files? Y/N [N]: ';
  deleteFiles = input(prompt, 's');
  if isempty(deleteFiles) str = 'N'; end;

  if strcmp(deleteFiles, 'Y')
    fprintf('Deleting %d individual geometry files (%d geometries contained)...\n', nFiles, nGeometries)
    for i = indices
      filename = strcat(filenamePrefix, '_G', sprintf('%05d',i), '.log');
      delete(filename);
    end
  end
end
