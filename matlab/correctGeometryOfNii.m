
% Start with the SIGMA template masked in-vivo image .
  S = load_untouch_nii('d0F12_S8_10m_UTE.nii');

% % load the t2w image, and the mask
%   t2w = load_untouch_nii('d12F21_S3_T2w_RARE_dcm_ic.nii');
%   mask = readnd('/mnt/data1/invicro/recon/2023-11-07_F21/d12F21_mask_S3_T2w.nd');

% % apply the mask, and orient the image volume so that it has the same orientation as the SIGMA image
%   t2w.img = flipdim( permute( t2w.img .* int16(mask), [ 1 3 2]), 3 );

% clone the SIGMA template, and replace the cloned image data with the t2w sampled image
  S2 = S;
  S2.img = t2w.img;

% now fix the header information...
  S2.hdr.dime.dim(2:4) = size(t2w.img);
  S2.hdr.dime.pixdim(2:4) = t2w.hdr.dime.pixdim(2:4);
  S2.hdr.hist.qoffset_x = -abs(t2w.hdr.hist.qoffset_x);
  S2.hdr.hist.qoffset_y = -abs(t2w.hdr.hist.qoffset_y);
  S2.hdr.hist.qoffset_z = -abs(t2w.hdr.hist.qoffset_z);
  S2.hdr.hist.srow_x = [ S2.hdr.dime.pixdim(2) 0 0 S2.hdr.hist.qoffset_x ];
  S2.hdr.hist.srow_y = [ 0 S2.hdr.dime.pixdim(3) 0 S2.hdr.hist.qoffset_y ];
  S2.hdr.hist.srow_z = [ 0 0 S2.hdr.dime.pixdim(4) S2.hdr.hist.qoffset_z ];
  S2.hdr.hist.descrip = 'T2w RARE image';

% save the t2w image using the new header information
  save_untouch_nii(S2,'t2w_test.nii');

% view in fsleyes, to confirm header is correct
  fsleyes  SIGMA_InVivo_Brain_Template_Masked.nii  t2w_test.nii

% adjust qoffset_? to put the origin point on the AC of the new t2w nii file


% repeat for the UTE data, using the same qoffset as the t2w (but different pixdim and datatype values)
  PC = S2;
  PC.img = flipdim(flipdim( permute( postc.img(:,:,:,1), [ 1 3 2 ]), 3 ),1);
  S2.hdr.dime.dim(2:4) = size(PC.img);
  S2.hdr.dime.pixdim(2:4) = postc.hdr.dime.pixdim(2:4);

  PC.hdr.dime.datatype = 16;
  PC.hdr.dime.bitpix = 32;
  PC.hdr.hist.descrip = '3D UTE post-contrast 10m40s vol-01';
  save_untouch_nii(PC,'pc_test.nii');
  
% run the registration
/usr/src/pkgs/ANTsX-ANTs-ae468a8/antsbin/ANTS-build/Examples//antsRegistration --verbose 1 --dimensionality 3 --float 0 --collapse-output-transforms 1 --output [ output,outputWarped.nii.gz,outputInverseWarped.nii.gz ] --interpolation Linear --use-histogram-matching 1 --winsorize-image-intensities [ 0.005,0.995 ] --initial-moving-transform [ t2w_test.nii,SIGMA_InVivo_Brain_Template_Masked.nii,1 ] --transform Rigid[ 0.1 ] --metric MI[ t2w_test.nii,SIGMA_InVivo_Brain_Template_Masked.nii,1,32,Regular,0.25 ] --convergence [ 200x100x50x20,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[ 0.1 ] --metric MI[ t2w_test.nii,SIGMA_InVivo_Brain_Template_Masked.nii,1,32,Regular,0.25 ] --convergence [ 200x100x50x20,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform SyN[ 0.1,3,0 ] --metric MI[ t2w_test.nii,SIGMA_InVivo_Brain_Template_Masked.nii,1,32] --convergence [ 200x140x100x50,1e-6,10 ] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox

% apply the estimated transforms to the labels

  antsApplyTransforms -i /code/analysis/SIGMA_Complete/SIGMA_Rat_Brain_Atlases/SIGMA_Functional_Atlas/SIGMA_Functional_Brain_Atlas_InVivo_Anatomical_Template.nii -n GenericLabel -o tmp.nii  -d 3 -v -t output1Warp.nii.gz -t output0GenericAffine.mat   -r t2w_test.nii