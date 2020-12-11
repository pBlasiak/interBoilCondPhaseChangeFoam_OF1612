# This file is NOT licensed under the GPLv3, which is the license for the rest
# of YouCompleteMe.
#
# Here's the license text for this file:
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>

import os
import ycm_core


foam_app = os.environ['FOAM_APP']
foam_src = os.environ['FOAM_SRC']
pb_proj = os.environ['WM_PROJECT_USER_DIR']
pb_sim = pb_proj + '/applications/solvers/heatTransfer/buoyantBoussinesqSuperFluidPimpleFoam/simplifiedSuperFluidModels/lnInclude'
pb_HeIIFoam = pb_proj + '/applications/solvers/solvers/HeIIFoam'
#drift_foam_dir = foam_app + '/multiphase/driftFluxFoam'
flags = [
    '-DSTF_EXPORTS',
    '-m64',
    '-Dlinux64',
    '-DWM_ARCH_OPTION=64',
    '-DWM_DP',
    '-DNoRepository',
    '-DWM_LABEL_SIZE=32', '-std=c++11',
    '-Wno-unused-local-typedefs',
    '-Wall',
    '-Wextra',
    '-Wno-unused-parameter',
    '-Wold-style-cast',
    '-Wnon-virtual-dtor',
    '-Wno-unused-variable',
    # '-DNDEBUG',
    # '-DUSE_CLANG_COMPLETER',
    # THIS IS IMPORTANT! Without a "-std=<something>" flag, clang won't know
    # which language to use when compiling headers. So it will guess. Badly. So
    # C++ headers will be compiled as C headers. You don't want that so ALWAYS
    # specify a "-std=<something>".
    # For a C project, you would set this to something like 'c99' instead of
    # 'c++11'.
    '-std=c++11',
    # ...and the same thing goes for the magic -x option which specifies the
    # language that the files to be compiled are written in. This is mostly
    # relevant for c++ headers.
    # For a C project, you would set this to 'c' instead of 'c++'.
    '-x',
    'c++',
    # This path will only work on OS X, but extra paths that don't exist are
    # not harmful
#    '-I' + drift_foam_dir + '/incompressibleTwoPhaseInteractingMixture',
    '-I' + foam_src + '/fvOptions/lnInclude',
    '-I' + foam_src + '/finiteVolume/lnInclude',
    '-I' + foam_src + '/meshTools/lnInclude',
    '-I' + foam_src + '/OpenFOAM/lnInclude',
    '-I' + foam_src + '/transportModels/incompressible/lnInclude',
    '-I' + foam_src + '/TurbulenceModels/incompressible/lnInclude',
    '-I' + foam_src + '/thermophysicalModels/radiation/lnInclude',
    '-I' + pb_proj + '/src',
    '-I' + pb_proj + '/applications',
    '-I' + pb_sim,
    '-I' + pb_HeIIFoam,
    '-I./src/',
    '-I.',
    '-I' + foam_src + '/OSspecific/POSIX/lnInclude'
]

# Set this to the absolute path to the folder (NOT the file!) containing the
# compile_commands.json file to use that instead of 'flags'. See here for
# more details: http://clang.llvm.org/docs/JSONCompilationDatabase.html
#
# Most projects will NOT need to set this to anything; you can just change the
# 'flags' list of compilation flags. Notice that YCM itself uses that approach.
compilation_database_folder = ''

if compilation_database_folder:
    database = ycm_core.CompilationDatabase(compilation_database_folder)
else:
    database = None


def DirectoryOfThisScript():
    return os.path.dirname(os.path.abspath(__file__))


def MakeRelativePathsInFlagsAbsolute(flags, working_directory):
    if not working_directory:
        return list(flags)
    new_flags = []
    make_next_absolute = False
    path_flags = ['-isystem', '-I', '-iquote', '--sysroot=']
    for flag in flags:
        new_flag = flag

        if make_next_absolute:
            make_next_absolute = False
            if not flag.startswith('/'):
                new_flag = os.path.join(working_directory, flag)

        for path_flag in path_flags:
            if flag == path_flag:
                make_next_absolute = True
                break

            if flag.startswith(path_flag):
                path = flag[len(path_flag):]
                new_flag = path_flag + os.path.join(working_directory, path)
                break

        if new_flag:
            new_flags.append(new_flag)
    return new_flags


def FlagsForFile(filename):
    if database:
        # Bear in mind that compilation_info.compiler_flags_ does NOT return a
        # python list, but a "list-like" StringVec object
        compilation_info = database.GetCompilationInfoForFile( filename )
        final_flags = MakeRelativePathsInFlagsAbsolute(
            compilation_info.compiler_flags_,
            compilation_info.compiler_working_dir_ )

    else:
        relative_to = DirectoryOfThisScript()
        final_flags = MakeRelativePathsInFlagsAbsolute( flags, relative_to )

    return {
        'flags': final_flags,
        'do_cache': True}

