
#include "link_nodes.h"

#include "system_structures.h"
#include "system.h"

#include "functor_de_link_nodes.h"
#include "functor_link_nodes.h"


void link_nodes(
	NodeInfoVecs& nodeInfoVecs,
	EdgeInfoVecs& edgeInfoVecs,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

		//Default fill values at 0.
		thrust::fill(nodeInfoVecs.links_made_individual_thread.begin(),
			nodeInfoVecs.links_made_individual_thread.end(), 0);

		thrust::fill(nodeInfoVecs.id_temp_linked_left.begin(),
				nodeInfoVecs.id_temp_linked_left.end(), 0);

		thrust::fill(nodeInfoVecs.id_temp_linked_right.begin(),
				nodeInfoVecs.id_temp_linked_right.end(), 0);

		thrust::counting_iterator<unsigned> counter(0);

		thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						counter,
						auxVecs.id_bucket_net_intc.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						counter,
						auxVecs.id_bucket_net_intc.begin())) + generalParams.max_node_count,
				nodeInfoVecs.links_made_individual_thread.begin(),//output
			functor_link_nodes(
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_y.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_collagen.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_elastin.data()),

				thrust::raw_pointer_cast(edgeInfoVecs.current_node_edge_count_vec.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_neighbors.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_length_zero.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_isedge_collagen.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_isedge_elastin.data()),

				thrust::raw_pointer_cast(auxVecs.id_value_expanded_net_intc.data()),
				thrust::raw_pointer_cast(auxVecs.key_begin_net_intc.data()),
				thrust::raw_pointer_cast(auxVecs.key_end_net_intc.data()),

				edgeInfoVecs.collagen_diameter,
				edgeInfoVecs.elastin_diameter,
				generalParams.max_nbr_count,
				generalParams.max_node_count,

				generalParams.max_links_per_iteration,
				thrust::raw_pointer_cast(nodeInfoVecs.id_temp_linked_left.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.id_temp_linked_right.data()) ) );


			/*
			for (unsigned i = 0; i < nodeInfoVecs.id_temp_linked_left.size(); i++) {
				unsigned varL = nodeInfoVecs.id_temp_linked_left[i];
				unsigned varR = nodeInfoVecs.id_temp_linked_right[i];

				if ((varL != 0) || (varR != 0))
					std::cout<< varL << " " <<varR << std::endl;
			}
			for (unsigned i = 0; i < edgeInfoVecs.global_neighbors.size(); i++) {
				if ((i > 0) && (i % ( generalParams.max_nbr_count) == 0)){
					std::cout << " " << std::endl;
				}
				unsigned varL = edgeInfoVecs.global_neighbors[i];
				if (varL < generalParams.max_node_count){
					std::cout << " " << varL;
				}
			}
			std::cout << " " << std::endl;
			std::cout<< nodeInfoVecs.node_loc_x[0] << " " << nodeInfoVecs.node_loc_y[0] << " " << nodeInfoVecs.node_loc_z[0]<< std::endl;
			std::cout<< nodeInfoVecs.node_loc_x[1] << " " << nodeInfoVecs.node_loc_y[1] << " " << nodeInfoVecs.node_loc_z[1]<< std::endl;
			*/
		/*unsigned begin = 479 * generalParams.max_nbr_count;
		unsigned end = begin + generalParams.max_nbr_count;
		for (unsigned i = begin; i < end; i++){
			unsigned id = edgeInfoVecs.global_neighbors[i];
			if (id < generalParams.max_node_count){
				std::cout<<" 479: "<< id <<std::endl;
			}
		}
		begin = 1004 * generalParams.max_nbr_count;
		end = begin + generalParams.max_nbr_count;
		for (unsigned i = begin; i < end; i++){
			unsigned id = edgeInfoVecs.global_neighbors[i];
			if (id < generalParams.max_node_count){
				std::cout<<" 1004: "<< id <<std::endl;
			}
		}*/
	/*	thrust::counting_iterator<unsigned> counterDeLink(0);

		thrust::transform(
						counterDeLink,
						counterDeLink + generalParams.max_node_count,
				nodeInfoVecs.delinksThreadMade.begin(),
			functor_de_link_nodes(
				thrust::raw_pointer_cast(edgeInfoVecs.global_neighbors.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_length_zero.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.current_node_edge_count_vec.data()),
				generalParams.max_nbr_count,
				generalParams.max_node_count ) );
*/

		//add links made by above function. links are double counted since threads create links from a->b and b->a
		unsigned num_placed_links = ceil( thrust::reduce(nodeInfoVecs.links_made_individual_thread.begin(),
			nodeInfoVecs.links_made_individual_thread.end(), 0, thrust::plus<unsigned>()) / 2 );

		//std::cout << " num_placed_links: " << num_placed_links << std::endl;


		//sort by increasing. Notice, the sorting must take place for the entire vector since threads write to different places
		thrust::sort_by_key(
			nodeInfoVecs.id_temp_linked_left.begin(),nodeInfoVecs.id_temp_linked_left.end(),
			nodeInfoVecs.id_temp_linked_right.begin(), thrust::greater<unsigned>() );

		thrust::stable_sort_by_key(
			nodeInfoVecs.id_temp_linked_right.begin(),nodeInfoVecs.id_temp_linked_right.end(),
			nodeInfoVecs.id_temp_linked_left.begin(), thrust::greater<unsigned>() );


		//the copy is not needed for the full vector. The only portion copied is 2 * num_placed_links since that many id's were set.
		thrust::copy(nodeInfoVecs.id_temp_linked_right.begin(),
			nodeInfoVecs.id_temp_linked_right.begin() + 2 * num_placed_links,
			nodeInfoVecs.host_id_right.begin() );
		thrust::copy(nodeInfoVecs.id_temp_linked_left.begin(),
			nodeInfoVecs.id_temp_linked_left.begin() + 2 * num_placed_links,
			nodeInfoVecs.host_id_left.begin());



		//old code, keep here in case of issues. use for validation.
		unsigned idL_init = nodeInfoVecs.host_id_left[0];
		unsigned idR_init = nodeInfoVecs.host_id_right[0];


		unsigned count = 0;
		//std::cout << "nodeInfoVecs.id_temp_linked_left.size(): " << nodeInfoVecs.id_temp_linked_left.size() << std::endl;
//		for (unsigned i = 1; i < nodeInfoVecs.id_temp_linked_left.size(); i++) {
//			//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and
//			//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
//			unsigned idL = nodeInfoVecs.host_id_left[i];
//			unsigned idR = nodeInfoVecs.host_id_right[i];
//
//			if ((idL == idL_init) && (idR == idR_init)){
//				count +=1;
//			}
//			else {
//				count = 0;
//			}
//			//reset initial id's
//			idL_init = idL;
//			idR_init = idR;
//
//
//			if ( ((idL != 0) || (idR != 0) ) && (count == 1)) {
//
//				//count edges
//				std::cout<<"placing id: "<< idL<<" " << idR<<std::endl;
//
//
//				nodeInfoVecs.host_edge_left[generalParams.current_edge_count] = (idL);
//				nodeInfoVecs.host_edge_right[generalParams.current_edge_count] = (idR);
//				generalParams.current_edge_count += 1;
//			}
//
//		}
		//end old code

		idL_init = nodeInfoVecs.host_id_left[0];
		idR_init = nodeInfoVecs.host_id_right[0];
		count = 0;
		for (unsigned i = 1; i < 2 * num_placed_links; i++) {
			//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and
			//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
			unsigned idL = nodeInfoVecs.host_id_left[i];
			unsigned idR = nodeInfoVecs.host_id_right[i];

			if ((idL == idL_init) && (idR == idR_init)) {
				count += 1;
			}
			else {
				count = 0;
			}
			//reset initial id's
			idL_init = idL;
			idR_init = idR;


			if (((idL != 0) || (idR != 0)) && (count == 1)) {

				//count edges
				//std::cout << "placing id from tester: " << idL << " " << idR << std::endl;
				//std::cout << " total edge count " << generalParams.current_edge_count << std::endl;
				////std::cout<< nodeInfoVecs.node_loc_x[idL] << " " << nodeInfoVecs.node_loc_y[idL] << " " << nodeInfoVecs.node_loc_z[idL]<< std::endl;
				////std::cout<< nodeInfoVecs.node_loc_x[idR] << " " << nodeInfoVecs.node_loc_y[idR] << " " << nodeInfoVecs.node_loc_z[idR]<< std::endl;
				//double dist = sqrt(
				//	(nodeInfoVecs.node_loc_x[idL] - nodeInfoVecs.node_loc_x[idR]) * (nodeInfoVecs.node_loc_x[idL] - nodeInfoVecs.node_loc_x[idR])+
				//	(nodeInfoVecs.node_loc_y[idL] - nodeInfoVecs.node_loc_y[idR]) * (nodeInfoVecs.node_loc_y[idL] - nodeInfoVecs.node_loc_y[idR])+
				//	(nodeInfoVecs.node_loc_z[idL] - nodeInfoVecs.node_loc_z[idR]) * (nodeInfoVecs.node_loc_z[idL] - nodeInfoVecs.node_loc_z[idR]));
				//std::cout<< "distance: " << dist << std::endl;
				nodeInfoVecs.host_edge_left[generalParams.current_edge_count] = (idL);
				nodeInfoVecs.host_edge_right[generalParams.current_edge_count] = (idR);
				generalParams.current_edge_count += 1;
			}

		}

	/*	unsigned globalcount = thrust::count_if(edgeInfoVecs.global_neighbors.begin(), edgeInfoVecs.global_neighbors.end(), is_less_than(generalParams.max_node_count));

		unsigned linksmade = *(thrust::max_element(links_made_individual_thread.begin(), links_made_individual_thread.end() ));
		unsigned delinksmade = *(thrust::max_element(delinksThreadMade.begin(), delinksThreadMade.end() ));
		std::cout<<"max links made this iteration: "<< linksmade << std::endl;
		std::cout<<"max unlinks made this iteration: "<< delinksmade << std::endl;

		std::cout<<"current_edge_count var: "<< generalParams.current_edge_count << std::endl;
		std::cout<<"current_edge_count global "<< globalcount/2 << std::endl;

		unsigned temp= thrust::reduce(	edgeInfoVecs.current_node_edge_count_vec.begin(),
			edgeInfoVecs.current_node_edge_count_vec.end());
		std::cout<<"current_edge_count dev: "<< temp << std::endl;
	*/




};
