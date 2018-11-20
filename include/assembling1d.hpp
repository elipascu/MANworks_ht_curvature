/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   assembling1d.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Miscelleanous assembly routines for the 1D network problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_HPP_
#define M3D1D_ASSEMBLING_1D_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>

namespace getfem {

//! Build the mass and divergence matrices for the 1D Poiseuille's problem,
//! @f$ M = \int_{\Lambda} c~u~v~ds @f$ and
//! @f$ D = \int_{\Lambda} \nabla u \cdot \mathbf{\lambda}\,p~ds @f$
/*!
	@param M         Computed mass matrix
	@param D         Computed divergence matrix
	@param mim       The integration method to be used
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_p      The finite element method for the pressure @f$ p @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for M
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param rg        The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_network_poiseuille
	(MAT & M, MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_fem & mf_data,
	 const VEC & coef,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_p.get_qdim() == 1 && mf_u.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
	// Build the local mass matrix Mvvi
	getfem::asm_mass_matrix_param(M, mim, mf_u, mf_data, coef, rg);
	// Build the local divergence matrix Dvvi
	generic_assembly 
	assem("l1=data$1(#3); l2=data$2(#3); l3=data$3(#3);"
		  "t=comp(Base(#2).Grad(#1).Base(#3));"
		  "M$1(#2,#1)+=t(:,:,1,i).l1(i)+t(:,:,2,i).l2(i)+t(:,:,3,i).l3(i);");
	assem.push_mi(mim);
	assem.push_mf(mf_u);
	assem.push_mf(mf_p);
	assem.push_mf(mf_data);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);
	assem.push_mat(D);
	assem.assembly(rg);
}

//! Build the mass and divergence matrices for the 1D Poiseuille's problem for compliant vessels,
//! @f$ M = \int_{\Lambda} c~u~v~ds @f$ and
//! @f$ D = \int_{\Lambda} \nabla u \cdot \mathbf{\lambda}\,p~ds @f$
/*!
	@param M         Computed mass matrix
	@param D         Computed divergence matrix
	@param mim       The integration method to be used
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_p      The finite element method for the pressure @f$ p @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for M
	@param coef      The coefficient for D
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param rg        The region where to integrate

	@ingroup asm
 */

template<typename MAT, typename VEC>
void
asm_network_poiseuille_rvar
(MAT & M, MAT & D,
	const mesh_im & mim,
	const mesh_fem & mf_u,
	const mesh_fem & mf_p,
	const mesh_fem & mf_data,
	const VEC & coefM,
	const VEC & coefD,
	const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	const mesh_region & rg = mesh_region::all_convexes()
)
{
	GMM_ASSERT1(mf_p.get_qdim() == 1 && mf_u.get_qdim() == 1,
		"invalid data mesh fem (Qdim=1 required)");
	// Build the local mass matrix Mvvi
	getfem::asm_mass_matrix_param(M, mim, mf_u, mf_data, coefM, rg);
	// Build the local divergence matrix Dvvi
	generic_assembly
		assem("l1=data$1(#3); l2=data$2(#3); l3=data$3(#3); cs=data$4(#2);"
			"t=comp(Base(#2).Grad(#1).Base(#3).Base(#2));"
			"t2=comp(Base(#2).Base(#1).Base(#3).Grad(#2));"
			"M$1(#2,#1)+=t(:,:,1,i,j).l1(i).cs(j)+t(:,:,2,i,j).l2(i).cs(j)+t(:,:,3,i,j).l3(i).cs(j)+ t2(:,:,i,j,1).l1(i).cs(j)+t2(:,:,i,j,2).l2(i).cs(j)+t2(:,:,i,j,3).l3(i).cs(j);");


	// vector_type ones(mf_p.nb_dof(),coefD[coefD.size()]);
	vector_type areap1(mf_p.nb_dof());
	//vector_type ones_p0(mf_data.nb_dof(),0.008);
	getfem::interpolation(mf_data, mf_p, coefD, areap1,  0);
	//getfem::interpolation(mf_p, mf_data, ones, ones_p0,  0);

	assem.push_mi(mim);
	assem.push_mf(mf_u);
	assem.push_mf(mf_p);
	assem.push_mf(mf_data);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);

	//  assem.push_data(ones);
	assem.push_data(areap1);
    //assem.push_data(ones_p0);
	assem.push_mat(D);              // output matrix
	assem.assembly(rg);
	vtk_export exp("./vtk/area"+std::to_string(rg.id())+".vtk");
	exp.exporting(mf_p);
    exp.write_point_data(mf_p, areap1, "r"); // write a scalar field
	// vtk_export exp2("./vtk/radp0"+std::to_string(rg.id())+".vtk");
	// exp2.exporting(mf_data);
    // exp2.write_cell_data( ones_p0, "r"); // write a scalar field
    // exp.write_cell_data(ones_p0, "rp0"); // write a p0 field

/*
	mesh_region mr_internal_face = inner_faces_of_mesh(mf_data.linked_mesh(),rg);
	vector_type coefDnext(coefD.size());
	vector_type coefDprev(coefD.size());

	for ( size_type i=1; i < coefD.size(); i++){
		coefDprev[i] =-coefD[i-1] + coefD[i];
     	cout << " area_diff["<<i<<"]  =  "<< coefDnext[i] << endl;
	}
	coefDnext[coefD.size()-1]=coefDnext[coefD.size()-2];

	coefDprev[0]=coefDprev[1];
	
		generic_assembly 
		assem_nodes("a=data$1(#3);"
			"t=comp(Base(#2).Base(#1).Base(#3));"
			"M$1(#2,#1)+=t(:,:,i).a(i);");
	assem_nodes.push_mi(mim);
	assem_nodes.push_mf(mf_u);
	assem_nodes.push_mf(mf_p);
	assem_nodes.push_mf(mf_data);
	assem_nodes.push_data(coefDnext);
	assem_nodes.push_mat(D);
	assem_nodes.assembly(mr_internal_face);
	*/

/*
	generic_assembly 
		assem_nodes("a=data$1(#3); anext=data$2(#3);aprev=data$3(#3);"
			"t=comp(Base(#2).Base(#1).Base(#3));"
			"t2=comp(Base(#2).Grad(#1).Base(#3));"
			"M$1(#2,#1)+=t(:,:,i).anext(i)*0.+ t2(:,:,1,i).anext(i)*0;");
	assem_nodes.push_mi(mim);
	assem_nodes.push_mf(mf_u);
	assem_nodes.push_mf(mf_p);
	assem_nodes.push_mf(mf_data);
	assem_nodes.push_data(coefD);
	assem_nodes.push_data(coefDnext);
	assem_nodes.push_data(coefDprev);
	assem_nodes.push_mat(D);
	assem_nodes.assembly(mr_internal_face);
*/

}


/*! Build the mixed boundary conditions for Poiseuille's problem
    @f$ M=\int_{\mathcal{E}_u} \frac{1}{\beta}\,u\,v~d\sigma@f$ and
    @f$ F=-\int_{\mathcal{E}_u} p0\,v~d\sigma-\int_{\mathcal{E}_p} g\,v~d\sigma@f$
 */
/*!
	@param M        BC contribution to Poiseuille's mass matrix
	@param F        BC contribution to Poiseuille's rhs
	@param mim      The integration method to be used
	@param mf_u     The finite element method for the velocity @f$u@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param P0       Array of values of the external pressure @f$p_0@f$
	@param R        Network radii
	@param rg       The region where to integrate
	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_network_bc
	(MAT & M, VEC & F,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const VEC & P0,
         const VEC & R
         //,const scalar_type beta
	 ) 
{
	// Aux data
	std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);

	for (size_type bc=0; bc < BC.size(); bc++) {

		size_type i = abs(BC[bc].branches[0]);
		size_type start = i*mf_u[i].nb_dof();
		scalar_type Ri = compute_radius(mim, mf_data, R, i);

		if (BC[bc].label=="DIR") { // Dirichlet BC
			// Add gv contribution to Fv
			scalar_type Rbc = compute_radius(mim, mf_data, R, i);//R[mf_data.ind_basic_dof_of_element(BC[bc].idx)[0]];
			scalar_type BCVal = BC[bc].value*pi*Rbc*Rbc;  // valore al bordo * area
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(ones, BCVal), BC[bc].rg);
		} 
		else if (BC[bc].label=="MIX") { // Robin BC
			// Add correction to Mvv
			MAT Mi(mf_u[i].nb_dof(), mf_u[i].nb_dof());
			getfem::asm_mass_matrix(Mi,
				mim, mf_u[i], BC[bc].rg);
			if (BC[bc].value == 0 ) GMM_WARNING1("You wanted to divide by BC[bc].value = 0 in asm_network_bc ");
			// dead ends are set as MIX with value 0
                        gmm::scale(Mi, pi*pi*Ri*Ri*Ri*Ri/BC[bc].value);
			gmm::add(gmm::scaled(Mi, -1.0), 
				gmm::sub_matrix(M,
					gmm::sub_interval(start, mf_u[i].nb_dof()),
					gmm::sub_interval(start, mf_u[i].nb_dof())));
			gmm::clear(Mi);	
			// Add p0 contribution to Fv
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(P0, pi*Ri*Ri), BC[bc].rg);			
		}
		else if (BC[bc].label=="INT") { // Internal Node
			GMM_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			GMM_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
		}
	}

}

/*! Build the mixed boundary conditions for Poiseuille's problem for compliant vessels
    @f$ M=\int_{\mathcal{E}_u} \frac{1}{\beta}\,u\,v~d\sigma@f$ and
    @f$ F=-\int_{\mathcal{E}_u} p0\,v~d\sigma-\int_{\mathcal{E}_p} g\,v~d\sigma@f$
 */
/*!
	@param M        BC contribution to Poiseuille's mass matrix
	@param F        BC contribution to Poiseuille's rhs
	@param mim      The integration method to be used
	@param mf_u     The finite element method for the velocity @f$u@f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param P0       Array of values of the external pressure @f$p_0@f$
	@param R        Network areas
	@param rg       The region where to integrate
	@ingroup asm
 */ 
template<typename VEC>
void
asm_network_bc_rvar
	(VEC & F,
	const mesh_im & mim,
	const std::vector<mesh_fem> & mf_u,
	const mesh_fem & mf_data,
	const std::vector<getfem::node> & BC,
	const VEC & P0,
    const VEC & area
	) 
{
	// Aux data
	std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);

	for (size_type bc=0; bc < BC.size(); bc++) {

		size_type i = abs(BC[bc].branches[0]);
		size_type start = i*mf_u[i].nb_dof();
		scalar_type area_loc = 0;

		if (BC[bc].label=="DIR") { // Dirichlet BC
			// Add gv contribution to Fv
			size_type k = mf_data.linked_mesh().convex_to_point(BC[bc].idx)[0];
			area_loc = area[mf_data.ind_basic_dof_of_element(k)[0]]; // also this works only for P0 data on vessels
			//cout << " asm bc   area_loc = " << area_loc << ",    node index = " << BC[bc].idx << endl;
			scalar_type BCVal = BC[bc].value*area_loc;  // valore al bordo * area
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(ones, BCVal), BC[bc].rg);
		} 
		/* else if (BC[bc].label=="MIX") { // Robin BC
			// Add correction to Mvv
			MAT Mi(mf_u[i].nb_dof(), mf_u[i].nb_dof());
			getfem::asm_mass_matrix(Mi,
				mim, mf_u[i], BC[bc].rg);
                        gmm::scale(Mi, pi*pi*Ri*Ri*Ri*Ri/BC[bc].value);
			gmm::add(gmm::scaled(Mi, -1.0), 
				gmm::sub_matrix(M,
					gmm::sub_interval(start, mf_u[i].nb_dof()),
					gmm::sub_interval(start, mf_u[i].nb_dof())));
			gmm::clear(Mi);	
			// Add p0 contribution to Fv
			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(start,mf_u[i].nb_dof())), 
				mim, mf_u[i], mf_data, gmm::scaled(P0, pi*Ri*Ri), BC[bc].rg);			
		}*/
		else if (BC[bc].label=="INT") { // Internal Node
			GMM_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			GMM_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Compliant vessels are solved only for DIR conditions. Unknown Boundary Condition"<< BC[bc].label << endl);
		}
	}

}



/*!
	Compute the network junction matrix @f$J=\langle[u],p\rangle_{\Lambda}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void
asm_network_junctions
	(MAT & J,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_p,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & J_data,
	 const VEC & radius
	 ) 
{
	GMM_ASSERT1 (mf_p.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1 (mf_u[0].get_qdim() == 1, 
		"invalid data mesh fem for velocity (Qdim=1 required)");
	GMM_ASSERT1 (getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK(1,0)" &&
		getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK_DISCONTINUOUS(1,0)",
		"invalid data mesh fem for pressure (k>0 required)");
	
	for (size_type i=0; i<mf_u.size(); ++i){ /* branch loop */

		//scalar_type Ri = compute_radius(mim, mf_data, radius, i);

		for (size_type j=0; j<J_data.size(); ++j){
			// Identify pressure dof corresponding to junction node
			VEC psi(mf_p.nb_dof());
			asm_basis_function(psi, mim, mf_p, J_data[j].rg);
			size_type row = 0;
			bool found = false;
			while (!found && row<mf_p.nb_dof()){
				found = (1.0 - psi[row] < 1.0E-06);
				if (!found) row++;
			}
			GMM_ASSERT1 (row!=0 && found,  // No junction in first point
				"Error in assembling pressure basis function");
			std::vector<long signed int>::const_iterator bb = J_data[j].branches.begin();
			std::vector<long signed int>::const_iterator be = J_data[j].branches.end();
			size_type last_, first_;
			vector_type dof_enum;
			int fine=0;
			for (getfem::mr_visitor mrv(mf_u[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto b : mf_u[i].ind_basic_dof_of_element(mrv.cv()))
				{dof_enum.emplace_back(b);
				fine++;}			
			first_=dof_enum[0];
			last_=dof_enum[fine-1];
			dof_enum.clear();
			scalar_type R_loc=0;
			//cout << " size convex to point "<< mf_data.linked_mesh().convex_to_point(J_data[j].idx).size() << endl;
			// mf_data.linked_mesh().region(i).convex_to_point(J_data[j].idx)[0] non puo farlo 
			for (auto k : mf_data.linked_mesh().convex_to_point(J_data[j].idx)){
				if ( mf_data.linked_mesh().region(i).is_in(k) ) {
					R_loc = radius[mf_data.ind_basic_dof_of_element(k)[0]]; // also this works only for P0 data on vessels
					//cout << " entro nell'if con dof = " << mf_data.ind_basic_dof_of_element(k)[0] << endl;
				}
			}
			//cout << "R_loc =  " << R_loc << ",  Ri =  "<< Ri<< endl;
			// Outflow branch contribution
			if (std::find(bb, be, i) != be){
				J(row, i*mf_u[i].nb_dof()+last_) -= pi*R_loc*R_loc;//col to be generalized!
				//cout << " primo if   R_loc =  " << R_loc << ",  Ri =  "<< Ri<< endl;
			}
			// Inflow branch contribution
			if (i!=0 && std::find(bb, be, -i) != be){
				J(row, i*mf_u[i].nb_dof()+first_) += pi*R_loc*R_loc;	//col to be generalized!
				//cout << " secondo if    R_loc =  " << R_loc << ",  Ri =  "<< Ri<< endl;
			}
		}
	}

} /* end of asm_junctions */

template<typename MAT, typename VEC>
void
asm_network_junctions_rvar
(MAT & J,
	const mesh_im & mim,
	const std::vector<mesh_fem> & mf_u,
	const mesh_fem & mf_p,
	const mesh_fem & mf_data,
	const std::vector<getfem::node> & J_data,
	const VEC & area
)
{
	GMM_ASSERT1(mf_p.get_qdim() == 1,
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1(mf_u[0].get_qdim() == 1,
		"invalid data mesh fem for velocity (Qdim=1 required)");
	GMM_ASSERT1(getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK(1,0)" &&
		getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK_DISCONTINUOUS(1,0)",
		"invalid data mesh fem for pressure (k>0 required)");

	for (size_type i = 0; i < mf_u.size(); ++i) { /* branch loop */

		for (size_type j = 0; j < J_data.size(); ++j) {
			// Identify pressure dof corresponding to junction node
			VEC psi(mf_p.nb_dof());
			asm_basis_function(psi, mim, mf_p, J_data[j].rg);
			size_type row = 0;
			bool found = false;
			while (!found && row < mf_p.nb_dof()) {
				found = (1.0 - psi[row] < 1.0E-06);
				if (!found) row++;
			}
			GMM_ASSERT1(row != 0 && found,  // No junction in first point
				"Error in assembling pressure basis function");
			std::vector<long signed int>::const_iterator bb = J_data[j].branches.begin();
			std::vector<long signed int>::const_iterator be = J_data[j].branches.end();
			size_type last_, first_;
			vector_type dof_enum;
			int fine = 0;
			for (getfem::mr_visitor mrv(mf_u[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
				for (auto b : mf_u[i].ind_basic_dof_of_element(mrv.cv()))
				{
					dof_enum.emplace_back(b);
					fine++;
				}
			first_ = dof_enum[0];
			last_ = dof_enum[fine - 1];
			dof_enum.clear();
			scalar_type area_loc = 0;
			//cout << " size convex to point "<< mf_data.linked_mesh().convex_to_point(J_data[j].idx).size() << endl;
			// mf_data.linked_mesh().region(i).convex_to_point(J_data[j].idx)[0] non puo farlo 
			for (auto k : mf_data.linked_mesh().convex_to_point(J_data[j].idx)) {
				if (mf_data.linked_mesh().region(i).is_in(k)) {
					area_loc = area[mf_data.ind_basic_dof_of_element(k)[0]]; // also this works only for P0 data on vessels
					//cout << " entro nell'if con dof = " << mf_data.ind_basic_dof_of_element(k)[0] << endl;
				}
			}
			//cout << "R_loc =  " << R_loc << ",  Ri =  "<< Ri<< endl;
			// Outflow branch contribution
			if (std::find(bb, be, i) != be) {
				J(row, i*mf_u[i].nb_dof() + last_) -= area_loc;//col to be generalized!
				//cout << " primo if   area_loc =  " << area_loc << "   ramo " << i << endl;
			}
			// Inflow branch contribution
			// if std::find finds the element, it returns an iterator to the element
			// if it doesn't find the element, it returns the value of the last element of the array
			// so if i=0 then it can only be inflow. if i=be, it works because it returns the iterator, not the value
			if (i != 0 && std::find(bb, be, -i) != be) {
				J(row, i*mf_u[i].nb_dof() + first_) += area_loc;	//col to be generalized!
				//cout << " secondo if    area_loc =  " << area_loc << "   ramo " << i  << endl;
			}
		}
	}

} /* end of asm_junctions_rvar */

} /* end of namespace */

#endif
