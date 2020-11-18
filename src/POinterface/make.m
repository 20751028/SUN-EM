% Makefile for mexRaySrv

%AUTHOR:
%
% J.C. Smit
%
%UPDATE RECORD:
%
% (2011-03-11) - Started

% --- User specified options ---
% Installation paths
if ispc
    %CUDA_install_path = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v3.2';
    CUDA_install_path = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1';
    %OptiX_install_path = 'C:\Program Files\NVIDIA Corporation\OptiX SDK 2.1.0';
    OptiX_install_path = 'C:\ProgramData\NVIDIA Corporation\OptiX SDK 6.5.0';
elseif ismac
    %CUDA_install_path = '/Developer/NVIDIA/CUDA-6.5';
    %CUDA_install_path = '/Developer/NVIDIA/CUDA-6.5';
    CUDA_install_path = '/Developer/NVIDIA/CUDA-7.5';
    OptiX_install_path = '/Developer/OptiX';
else
    CUDA_install_path = '/usr/local/cuda';
    %OptiX_install_path = '/home/jcsmit/NVIDIA-OptiX-SDK-2.1.0-linux32';
    %OptiX_install_path = '/home/jcsmit/NVIDIA-OptiX-SDK-2.1.0-linux64';
    %OptiX_install_path = '/home/jcsmit/NVIDIA-OptiX-SDK-2.1.1-linux64';
    %OptiX_install_path = '/home/jcsmit/NVIDIA-OptiX-SDK-2.5.1-linux64';
    %OptiX_install_path = '/home/jcsmit/NVIDIA-OptiX-SDK-2.6.0-linux64';
    %OptiX_install_path = '/home/jcsmit/Documents/SDK/NVIDIA-OptiX-SDK-2.6.0-linux64';
%     OptiX_install_path = '/opt/optix/current';
    OptiX_install_path = '/shared/apps/nvidia/optix/NVIDIA-OptiX-SDK-3.9.1-linux64/lib64';
end

% --- Automatic options ---
% Paths setup
CUDA_include_path = fullfile(CUDA_install_path,'include');
OptiX_include_path = fullfile(OptiX_install_path,'include');
OptiX_internal_include_path = fullfile(OptiX_include_path,'internal');
OptiX_prime_include_path = fullfile(OptiX_include_path,'optix_prime');
OptiX_lib_path = fullfile(OptiX_install_path,'lib64');

% Initialise mex arguments
args = {};

% Source files
args{end+1} = sprintf('-I%s',CUDA_include_path);
args{end+1} = sprintf('-I%s',OptiX_include_path);
args{end+1} = sprintf('-I%s',OptiX_internal_include_path);
args{end+1} = sprintf('-I%s',OptiX_prime_include_path);
args{end+1} = 'mexRaytracer.cpp';

% Linking options
args{end+1} = sprintf('-L%s',OptiX_lib_path);
if ispc
    args{end+1} = '-loptix.6.5.0';
   % args{end+1} = '-loptixu.1';
     args{end+1} = '-loptix_prime.6.5.0';
else
    if ismac
        args{end+1} = '-loptix';
        args{end+1} = '-loptixu';
        %     args{end+1} = '-loptix_prime';
    else
        args{end+1} = '-loptix';
        args{end+1} = '-loptixu';
        %     args{end+1} = '-loptix_prime';
    end
    args{end+1} = sprintf('LDFLAGS=$LDFLAGS -Wl,-rpath,%s',OptiX_lib_path);
end

% Debug
 %args{end+1} = '-g';

% Call mex with all arguments
mex(args{:});

% EOF
