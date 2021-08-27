# electrons density
A=0
for i in `seq 0.5 0.5 4`; do
    echo "electrons density ${A} ${i}" >> list
    ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'density\' time_range=[${A}e-8,${i}e-8] clim=[1e16,1.5e17] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
    A=${i}
done

# ions density 1
A=0
for i in `seq 0.5 0.5 2`; do
    echo "ions density ${A} ${i}" >> list
    ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'density\' time_range=[${A}e-8,${i}e-8] clim=[5e16,1.3e17] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
    A=${i}
done

# ions density 2
A=2
for i in `seq 2.5 0.5 4`; do
    echo "ions density ${A} ${i}" >> list
	./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'density\' time_range=[${A}e-8,${i}e-8] clim=[3e16,1.5e17] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
	A=${i}
    done
done

# electrons temperature
## special
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[0,5e-9] clim=[0,1e2] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[5e-9,1e-8] clim=[0,1e3] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[1e-8,1.5e-8] clim=[0,1e3] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[1.5e-8,2e-8] clim=[0,5e3] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[2e-8,2.5e-8] clim=[0,1e4] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'

## regular
A=2.5
for i in `seq 3 0.5 4`; do
    echo "electrons temperature ${A} ${i}" >> list
    ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[${A}e-8,${i}e-8] clim=[0,1e4] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
    A=${i}
done

# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[3e-8,3.5e-8] clim=[0,1e4] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'electrons\' cmp=\'temperature\' time_range=[3.5e-8,4e-8] clim=[0,1e4] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'

# temperature ions
## special
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[0,5e-9] clim=[0,0.5] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[5e-9,1e-8] clim=[0,1] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[1e-8,1.5e-8] clim=[0,5] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'

## regular
A=1.5
for i in `seq 2 0.5 4`; do
    echo "ions temperature ${A} ${i}" >> list
    ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[${A}e-8,${i}e-8] clim=[0,50] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
    A=${i}
done

# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[2e-8,2.5e-8] clim=[0,50] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[2.5e-8,3e-8] clim=[0,50] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[3e-8,3.5e-8] clim=[0,50] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch tools/base/mesh_movie.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/\' cmap=\'gray_r\' specie=\'ions\' cmp=\'temperature\' time_range=[3.5e-8,4e-8] clim=[0,50] frame_size=[0,0,300,1980] frame_step=2 image_interpolation=\'bilinear\'

# 3E beam density
## special
./tools/nbrun.sh ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-1.5e6,1.5e6] clim_z=[-1.5e6,1.5e6] clim_beam=[0,1e16] time_range=[0,0.2e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-3e6,3e6] clim_z=[-3e6,3e6] clim_beam=[0,5e15] time_range=[0.2e-8,0.4e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-4e6,4e6] clim_z=[-4e6,4e6] clim_beam=[0,5e15] time_range=[0.4e-8,1e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'

## regular 1
A=1
for i in `seq 1.5 0.5 2`; do
    echo "ions temperature ${A} ${i}" >> list
    ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-3e6,3e6] clim_z=[-3e6,3e6] clim_beam=[0,5e15] time_range=[${A}e-8,${i}e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
    # ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-3e6,3e6] clim_z=[-3e6,3e6] clim_beam=[0,5e15] time_range=[1.5e-8,2e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
    A=${i}
done


## regular 2
A=2
for i in `seq 2.5 0.5 4`; do
    echo "ions temperature ${A} ${i}" >> list
    ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-2e6,2e6] clim_z=[-2e6,2e6] clim_beam=[0,5e15] time_range=[${A}e-8,${i}e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'

# ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-2e6,2e6] clim_z=[-2e6,2e6] clim_beam=[0,5e15] time_range=[2e-8,2.5e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-2e6,2e6] clim_z=[-2e6,2e6] clim_beam=[0,5e15] time_range=[2.5e-8,3e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-2e6,2e6] clim_z=[-2e6,2e6] clim_beam=[0,5e15] time_range=[3e-8,3.5e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
# ./tools/nbrun.sh -batch ./tools/movie_3_E_beam_density.ipynb data_path=\'/media/nfs/PiCoPiC/m50s1/data.h5\' frame_size=[0,0,300,1980] frame_step=2 clim_r=[-2e6,2e6] clim_z=[-2e6,2e6] clim_beam=[0,5e15] time_range=[3.5e-8,4e-8] cmap=\'gray_r\' image_interpolation=\'bilinear\'
    A=${i}
done
