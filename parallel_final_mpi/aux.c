MPI_Datatype mpi_node_type;

int blocklengths[5] = {1,n_dims,1,1,1};
MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_LONG,MPI_LONG,MPI_LONG};
MPI_Aint     offsets[5];
offsets[0] = offsetof(node, radius);
offsets[1] = offsetof(node, center);
offsets[2] = offsetof(node, L);;
offsets[3] = offsetof(node, R);
offsets[4] = offsetof(node, id);
MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_node_type);
MPI_Type_commit(&mpi_node_type);