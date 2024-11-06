% update the nifti header description to include the acquisition time-stamp, the temporal resolution, and the subject, day, and image count

% to be run in a subdirectory below 'analysis', after the prepped
% nifti files have been burst into separate images


if ~exist('l','var')
    curdir = pwd;
    tmp = regexp( curdir, 'analysis');
    l = dir( [ curdir(1:tmp-1) 'prep' ] );
end;

% if ~exist('f','var');
f = dir('*.nii');
% end;

if ~exist('dbg','var'); dbg = 1; end;

for cnt=1:length(f)

    tmp = regexp(f(cnt).name,'(d\d+)([MF]\d+)_(S\d+)_(\d+m)_(.+?).nii','tokens');
    dd = ''; subj = ''; scan = ''; tres = ''; vnum = '';
    if ~isempty(tmp)
        dd = tmp{1}{1}; subj = tmp{1}{2}; scan = tmp{1}{3}; tres = tmp{1}{4};
        name = tmp{1}{5};
        tmp = regexp(name,'\d+$');
        if ~isempty(tmp)
            vnum = name(tmp:end);
            name = name(1:tmp-1);
        end;
    end
    indx = []; rpt = [];
    for  lcnt=1:length(l);
        if ( strfind( l(lcnt).name, subj ) );
            indx = [ indx lcnt] ;
            rptname = [ l(lcnt).folder filesep l(lcnt).name filesep dd subj '_scan_report.txt' ];
            if isfile( rptname )
                rpt = load_scan_rpt( rptname );
            end
        end;
        if ~isempty(rpt)
            break;
        end;
    end;

    if ~isempty(rpt)
        eval([ 'base_ts = rpt.' scan '.Time_Stamp;' ]);
        base_ts = datetime(base_ts);
        base_ts.Format = 'uuuu-MM-dd HH:mm:ss.SS';
        
        trfact = 1;
        if strfind( tres, '2m' )
            trfact = 4;
        elseif strfind( tres, '1m' )
            trfact = 8;
        end

        fprintf('%s  %s\n',f(cnt).name,vnum);
        majorloop = 0; minorloop = 0;
        if ~isempty(vnum)
            vnum = str2num(vnum);
            minorloop = mod(( vnum - 1 )/trfact,1);
            majorloop = (vnum - 1)/trfact - minorloop;
            % images without numbers are presumed to be pre-contrast
            cntrst = ' post-contrast';
        else
            % images without numbers are presumed to be pre-contrast
            cntrst = ' pre-contrast';
        end;

        eval([ 'minorenc = floor( rpt.' scan '.N_Encodes / trfact );' ]);
        eval([ 'nencs = rpt.' scan '.N_Encodes * majorloop + minorenc * minorloop * trfact; ']);
        eval([ 'esec  = nencs * rpt.' scan '.TR / 1000;' ]);
        if (dbg)
            fprintf('\tbase timestamp: %s\n',base_ts);
            fprintf('\tmajor loop: %d, minor loop: %d\n', majorloop, minorloop * trfact);
            fprintf('\telapsed sec:    %f\n', esec );
            fprintf('\tacq  timestamp: %s\n', datetime(base_ts) + seconds(esec) );
        end
        pc = load_untouch_nii( f(cnt).name );
        curdesc = pc.hdr.hist.descrip;
        tmp = findstr( curdesc, '::' );
        if isempty(tmp)
            if ( findstr( curdesc, subj ) )
                curdesc = sprintf('%s %s',string(datetime(base_ts) + seconds(esec) ),[ ':: ' curdesc ]);
            else
                eval([' mn = seconds( minorenc * rpt.' scan '.TR/1000 );']);
                res = sprintf(' %dm%ds ',floor(minutes(mn)), floor(seconds(mn) - 60*floor(minutes(mn))) );
                curdesc = sprintf('%s %s%d',string(datetime(base_ts) + seconds(esec) ),[ ':: ' subj ' ' dd res name cntrst ' image '], vnum  );
            end
        else
            curdesc = sprintf('%s %s',string(datetime(base_ts) + seconds(esec) ), curdesc(tmp:end) );
        end;
        if strfind(pc.hdr.hist.descrip,' realigned');
            curdesc = [ curdesc ' - realigned' ];
        end
        fprintf('  current descrip: %s\n',pc.hdr.hist.descrip);
        fprintf('      new descrip: %s\n',curdesc);

        if (dbg == 0)
            pc.hdr.hist.descrip = curdesc;
            save_untouch_nii( pc, pc.fileprefix );
        end;
        
    else
        fprintf([' scan report not found for ' f(cnt).name '\n']);
    end;
    
end
