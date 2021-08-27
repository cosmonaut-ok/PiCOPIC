
SPECIE=electrons
CMP=temperature

if [ ! -z ${1} ]; then
    CMP=$1
fi

if [ ! -z ${2} ]; then
    SPECIE=$2
fi

mpv tools/base/mesh_movie_${CMP}-${SPECIE}_0-5e-09.avi tools/base/mesh_movie_${CMP}-${SPECIE}_5e-09-1e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_1e-08-1.5e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_1.5e-08-2e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_2e-08-2.5e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_2.5e-08-3e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_3e-08-3.5e-08.avi tools/base/mesh_movie_${CMP}-${SPECIE}_3.5e-08-4e-08.avi
