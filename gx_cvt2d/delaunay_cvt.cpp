/*
*  _____ _____ ________  _
* /  __//  __//  __/\  \//
* | |  _|  \  |  \   \  / 
* | |_//|  /_ |  /_  /  \ 
* \____\\____\\____\/__/\\
*
* Graphics Environment for EXperimentations.
*  Copyright (C) 2006 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: 
*
*     ALICE Project - INRIA
*     INRIA Lorraine, 
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX 
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs. 
*/

#include "delaunay_cvt.h"
#include "geometry.h"
#include "lloyd_energy.h"
#include <glut_viewer/glut_viewer.h>

#include "generated/L.h"

#include "generated/L2.h"
#include "generated/L4.h"
#include "generated/L6.h"
#include "generated/L8.h"
#include "generated/L10.h"
#include "generated/L12.h"
#include "generated/L14.h"
#include "generated/L16.h"
#include "generated/L18.h"
#include "generated/L20.h"

#include "generated/Ltheta2.h"
#include "generated/Ltheta4.h"
#include "generated/Ltheta6.h"
#include "generated/Ltheta8.h"
#include "generated/Ltheta10.h"
#include "generated/Ltheta12.h"
#include "generated/Ltheta14.h"
#include "generated/Ltheta16.h"
#include "generated/Ltheta18.h"
#include "generated/Ltheta20.h"

#include "generated/Ltheta_rho2.h"
#include "generated/Ltheta_rho4.h"
#include "generated/Ltheta_rho6.h"
#include "generated/Ltheta_rho8.h"
#include "generated/Ltheta_rho10.h"
#include "generated/Ltheta_rho12.h"
#include "generated/Ltheta_rho14.h"
#include "generated/Ltheta_rho16.h"
#include "generated/Ltheta_rho18.h"
#include "generated/Ltheta_rho20.h"

#include <fstream>

static bool lbfgs_redraw = true ;

void funcgrad_cvt2d(int N, double* x, double& f, double* g);
void newiteration_cvt2d(int N, const double* x, double f, const double* g, double gnorm);


namespace Geex {

	DelaunayCVT* DelaunayCVT::instance_ = nil ;

	DelaunayCVT::DelaunayCVT(Delaunay* delaunay) : delaunay_(delaunay) {
		gx_assert(instance_ == nil) ;
		instance_ = this ;
		symbolic_ = false ;

		//      funcs_table_[0] = new GenericLloydFuncs<L> ;

		funcs_table_[0] = new GenericLloydFuncs<L2> ;
		funcs_table_[1] = new GenericLloydFuncs<L4> ;
		funcs_table_[2] = new GenericLloydFuncs<L6> ;
		funcs_table_[3] = new GenericLloydFuncs<L8> ;
		funcs_table_[4] = new GenericLloydFuncs<L10> ;
		funcs_table_[5] = new GenericLloydFuncs<L12> ;
		funcs_table_[6] = new GenericLloydFuncs<L14> ;
		funcs_table_[7] = new GenericLloydFuncs<L16> ;
		funcs_table_[8] = new GenericLloydFuncs<L18> ;
		funcs_table_[9] = new GenericLloydFuncs<L20> ;

		funcs_table_theta_[0] = new GenericLloydFuncs<Ltheta_rho2>(2) ;
		funcs_table_theta_[1] = new GenericLloydFuncs<Ltheta_rho4>(2) ;
		funcs_table_theta_[2] = new GenericLloydFuncs<Ltheta_rho6>(2) ;
		funcs_table_theta_[3] = new GenericLloydFuncs<Ltheta_rho8>(2) ;
		funcs_table_theta_[4] = new GenericLloydFuncs<Ltheta_rho10>(2) ;
		funcs_table_theta_[5] = new GenericLloydFuncs<Ltheta_rho12>(2) ;
		funcs_table_theta_[6] = new GenericLloydFuncs<Ltheta_rho14>(2) ;
		funcs_table_theta_[7] = new GenericLloydFuncs<Ltheta_rho16>(2) ;
		funcs_table_theta_[8] = new GenericLloydFuncs<Ltheta_rho18>(2) ;
		funcs_table_theta_[9] = new GenericLloydFuncs<Ltheta_rho20>(2) ;


		funcs_ = funcs_table_[0] ; 

		Lp_ = 0.0 ; 
		X_scale_ = 0.0 ;
		Y_scale_ = 0.0 ;

		use_theta_ = GL_FALSE ;
		mode_ = LLOYD ;
		aniso_mode_ = CONSTANT ;

		center_mode_ = CENTROID;

		has_affinity_ = false ;
		affinity_weight_ = 0.0 ;

		//dxy add: init
		direction_field_ = GL_TRUE;
		direction_energy_mode_ = ANGLE_DUAL;
		direction_field_mode_ = CONSTANT;
		direction_constrained_edge_ = GL_FALSE;
		direction_boundary_edge_ = GL_FALSE;
		direction_gradient_ = GL_FALSE;
		use_quad_field_ = GL_FALSE;
		use_true_grad_ = GL_FALSE;
		use_lloyd_energy_ = GL_TRUE;
		use_topo_optimization_ = GL_FALSE;

		direction_angle_ = 0.0;
		direction_factor_ = 1.0 /*0.0015*/;
		quad_factor_ = 0.5;
		///
	}

	DelaunayCVT::~DelaunayCVT() {
		for(int i=0; i<MAX_P+1; i++){
			delete funcs_table_[i];
			delete funcs_table_theta_[i];
		}
		instance_ = nil ; 
	}

	void DelaunayCVT::lloyd(int nb_iter, bool redraw) {
		for(unsigned int k=0; k<nb_iter; k++) {
			std::vector<vec2> new_points ;
			FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
				double V ;
				vec2 g ;
				if (center_mode_ == CENTROID)
					get_cell_centroid(v, g, V) ;
				else if(center_mode_ == QUASI_INCENTER)
					get_cell_quasi_incenter(v, g);
				new_points.push_back(g) ;
			}
			delaunay_->clear() ;
			delaunay_->begin_insert() ;
			for(unsigned int j=0; j<new_points.size(); j++) {
				delaunay_->insert(new_points[j]) ;
			}
			delaunay_->end_insert(redraw) ;
			double F = lloyd_energy() ;
			if(redraw) {
				std::cerr << "Lloyd energy = " << F << std::endl ;
			}
		}
	}


	//----------------------------------------------------------------------------------------------------

	void DelaunayCVT::newton_lloyd(int nb_iter, bool redraw) {

		use_theta_ = (mode() == THETA) ;

		if(use_theta_) {
			FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
				if(false && it->dual_intersects_boundary) {
					it->theta = default_theta(it) ;
				}
				it->rho = default_rho(it) ;
			}
		} 

		lbfgs_redraw = redraw ;

		int iLp = int(Lp_) ;
		iLp = gx_min(iLp, int(MAX_P)) ;
		iLp = gx_max(iLp, 0) ;

		if(use_theta_) {
			symbolic_ = true ;
			funcs_ = funcs_table_theta_[iLp] ;
		} else {
			if(iLp == 0) { 
				symbolic_ = false ; 
			} else {
				symbolic_ = true ;
				iLp-- ;
				funcs_ = funcs_table_[iLp] ;
			}
		}

		int n = use_theta_ ? delaunay_->nb_vertices() * 4 : delaunay_->nb_vertices() * 2 ;
		int m = 7 ;

		double* x = new double[n];


		if(use_theta_){
			for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
				Delaunay::Vertex_handle it = delaunay_->all_vertices_[i] ;
				x[4*i  ] = it->point().x() ;
				x[4*i+1] = it->point().y() ;
				x[4*i+2] = it->theta ;
				x[4*i+3] = it->rho ;
			}
		} else {
			for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
				Delaunay::Vertex_handle it = delaunay_->all_vertices_[i] ;
				x[2*i  ] = it->point().x() ;
				x[2*i+1] = it->point().y() ;
			}
		}

		double epsg = 0, epsf=0, epsx=0;

		//Optimizer* opt = new LBFGSOptimizer();
		Optimizer* opt = new HLBFGSOptimizer();

		opt->set_epsg(epsg);
		opt->set_epsf(epsf);
		opt->set_epsx(epsx);

		opt->set_M(m);
		opt->set_N(n);
		opt->set_max_iter(nb_iter);

		opt->set_newiteration_callback(newiteration_cvt2d);
		opt->set_funcgrad_callback(funcgrad_cvt2d);

		opt->optimize(x) ;

		set_vertices(x) ;
		delete opt;
		delete [] x;
	}

	void DelaunayCVT::set_vertices(const double* x) {
		std::vector<double> energy ;
		std::vector<bool> locked ;
		///dxy add
		std::vector<vec2> lloyd_grad ;
		std::vector<vec2> direction_grad ;
		std::vector<double> area ;
		///

		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			energy.push_back(it->energy) ;
			locked.push_back(it->locked) ;
			///dxy add
			lloyd_grad.push_back(it->lloyd_grad) ;
			direction_grad.push_back(it->direction_grad);
			area.push_back(it->area);
			///
		}

		int n = use_theta_ ? delaunay_->nb_vertices() * 4 : delaunay_->nb_vertices() * 2 ;

		delaunay_->clear() ;
		delaunay_->begin_insert() ;
		{
			int i = 0 ;
			unsigned int e_i = 0 ;
			while(i+1 < n) {
				Delaunay::Vertex_handle v = delaunay_->insert(vec2(x[i], x[i+1])) ;
				if(v == nil || e_i > energy.size()-1 || e_i > locked.size()-1) {
				} else {
					v->energy = energy[e_i] ;
					v->locked = locked[e_i] ;
					///dxy add
					v->lloyd_grad = lloyd_grad[e_i];
					v->direction_grad = direction_grad[e_i];
					v->area = area[e_i];
					///
				}
				if(v != nil && use_theta_) {
					v->theta = x[i+2] ;
					v->rho = x[i+3] ;
				}
				i+=2 ;
				if(use_theta_) {
					i+=2 ;
				} else {
					v->theta = default_theta(v) ;
				}
				e_i++ ;
			}
		}
		delaunay_->end_insert(false) ;
	}


	//----------------------------------------------------------------------------------------------------

	void DelaunayCVT::get_cell_centroid(Delaunay::Vertex_handle v, vec2& g, double& V) {
		vec2 p0 = to_geex(v->point()) ;
		g.x = 0.0 ; g.y = 0.0 ; V = 0.0 ;
		if(delaunay_->dual_cell_intersects_boundary(v)) { 
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			for(unsigned int i=0; i<P->size(); i++) {
				const vec2& p1 = (*P)[i].vertex[0] ;
				const vec2& p2 = (*P)[i].vertex[1] ;
				double Vi = triangle_area(p0, p1, p2) ;
				vec2 Gi = triangle_centroid(p0, p1, p2) ;
				V += Vi ;
				g += Vi * Gi ;
			}
		} else {
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			Delaunay::Face_circulator jt = it ; jt++ ;
			do {
				const vec2& p1 = it->dual ;
				const vec2& p2 = jt->dual ;
				double Vi = triangle_area(p0, p1, p2) ;
				vec2 Gi = triangle_centroid(p0, p1, p2) ;
				V += Vi ;
				g += Vi * Gi ;
				it++ ; jt++ ;
			} while(it != delaunay_->incident_faces(v)) ;
		}
		double s = (1.0 / V) ;
		g.x *= s ;
		g.y *= s ;
	}
	//----------------------------------------------------------------------------------------------------
	void DelaunayCVT::get_cell_quasi_incenter(Delaunay::Vertex_handle v, vec2& g) {
		vec2 p0 = to_geex(v->point()) ;
		double mat[2][2] ={{0,0},{0,0}}, b[2] = {0,0};  
		vec2 edge_dir, normal_dir;
		double tmp;
		if(delaunay_->dual_cell_intersects_boundary(v)) { 
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			for(unsigned int i=0; i<P->size(); i++) {
				const vec2& p1 = (*P)[i].vertex[0] ;
				const vec2& p2 = (*P)[i].vertex[1] ;
				edge_dir = p1 - p2;
				double len = edge_dir.length();
				if (len > 0)
				{
					normal_dir.x = edge_dir.y / len; normal_dir.y = -edge_dir.x / len;
					mat[0][0] += len * normal_dir.x * normal_dir.x;
					mat[0][1] += len * normal_dir.x * normal_dir.y;
					mat[1][1] += len * normal_dir.y * normal_dir.y;
					tmp = dot(p1 + p2, normal_dir) * 0.5;
					b[0] += len * tmp * normal_dir.x;
					b[1] += len * tmp * normal_dir.y;
				}
			}
		} else {
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			Delaunay::Face_circulator jt = it ; jt++ ;
			do {
				const vec2& p1 = it->dual ;
				const vec2& p2 = jt->dual ;
				edge_dir = p1 - p2;
				double len = edge_dir.length();
				if (len > 0)
				{
					normal_dir.x = edge_dir.y / len; normal_dir.y = -edge_dir.x / len;
					mat[0][0] += len * normal_dir.x * normal_dir.x;
					mat[0][1] += len * normal_dir.x * normal_dir.y;
					mat[1][1] += len * normal_dir.y * normal_dir.y;
					tmp = dot(p1 + p2, normal_dir) * 0.5;
					b[0] += len * tmp * normal_dir.x;
					b[1] += len * tmp * normal_dir.y;
				}
				it++ ; jt++ ;
			} while(it != delaunay_->incident_faces(v)) ;
		}

		mat[1][0] = mat[0][1];
		double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		if(det != 0)
		{
			g.x = (b[0]*mat[1][1] - b[1]*mat[0][1])/det;
			g.y = (b[1]*mat[0][0] - b[0]*mat[1][0])/det;
		}
		else
			g = p0;
	}

	///dxy add
	//get field angle and gradient at position v
	void DelaunayCVT::get_field(const vec2& v, double& angle, vec2& grad) {
		switch (direction_field_mode_)
		{
		case F_CONSTANT:
			//constant field
			{
				angle = direction_angle_; 
				grad = vec2(0.0, 0.0);
			}

			break;
		case F_CIRCLE:
			//circle field center(0.5, 0.5)
			{
				vec2 center(sqrt(0.5 * delaunay_->nb_vertices()), sqrt(0.5 * delaunay_->nb_vertices()));
				vec2 d = v - center;
				if (d.length() == 0)
				{
					angle = direction_angle_;
					grad = vec2(0.0, 0.0);
				}
				else {
					angle = acos(d.x/d.length());
					if (d.y < 0) {
						angle = 2 * M_PI - angle;
					}
					angle += 0.5 * M_PI + direction_angle_;	
					grad = vec2(-d.y/d.length2(), d.x/d.length2());
				}
			}
			break;
		case F_SIN:
			//(sin y, sin x) center(0.5, 0.5)
			{
				vec2 center(sqrt(0.5 * delaunay_->nb_vertices()), sqrt(0.5 * delaunay_->nb_vertices()));
				vec2 d = v - center;
				d /= sqrt(0.5 * delaunay_->nb_vertices());
				d *= M_PI;
				if (d.length() == 0)
				{
					angle = direction_angle_;
					grad = vec2(0.0, 0.0); //singularity

				}
				else {
					vec2 fvec(sin(d.y), sin(d.x));
					angle = acos(fvec.x / fvec.length());
					if (fvec.y < 0)
					{
						angle = 2 * M_PI - angle;
					}
					angle += direction_angle_;
					double s2 = 1.0 / (sin(d.y)*sin(d.y) + sin(d.x)*sin(d.x));
					grad = s2 * vec2(cos(d.x) * sin(d.y), - sin(d.x) * cos(d.y)); 
				}
			}
			break;
		default:
			break;
		}
	}
	///
	//----------------------------------------------------------------------------------------------------
	void DelaunayCVT::get_edge_direction_fg(const vec2& edge, const vec2& v, double& f, vec2& grad_f) {
		//assert: edge.length() > 0
		double beta = acos(edge.x / edge.length());
		if (edge.y < 0)
		{
			beta = 2 * M_PI - beta;
		}
		double field_angle;
		vec2 field_grad;
		get_field(v, field_angle, field_grad);
		if (use_true_grad_)
		{
			double field_angle1, field_angle2;
			vec2 field_grad1;
			get_field(v+edge, field_angle2, field_grad1);
			get_field(v, field_angle1, field_grad1);
			f = 0.5 * (1 - cos(6 * (beta - field_angle1)));
			vec2 grad_beta = vec2(edge.y, - edge.x);
			//grad_beta /= edge.length();
			grad_beta /= edge.length2();
			if (direction_gradient_)
			{
				grad_f = 3 * sin(6*(beta-field_angle1)) * (grad_beta-field_grad1) + 3 * sin(6*(beta-field_angle2)) * grad_beta;
			}
			else {
				grad_f = 3 * (sin(6*(beta-field_angle)) + sin(6*(beta-field_angle2))) * grad_beta;
			}
			return;
		}

		if (!use_quad_field_)
		{
			f = 0.5 * (1 - cos(6 * (beta - field_angle)));
			grad_f = vec2(edge.y, - edge.x);
			//grad_f /= edge.length();
			grad_f /= edge.length2();
			if (direction_gradient_)
			{
				grad_f -= field_grad;
			}
			grad_f = (3 * sin(6 * (beta-field_angle))) * grad_f;
		}
		else {
			double lambda = quad_factor_;
			double angle_diff  = beta - field_angle;
			f = 0.5 * (1 - lambda * cos(8 * angle_diff) + (lambda-1) * cos(4 * angle_diff));
			grad_f = vec2(edge.y, - edge.x);
			//grad_f /= edge.length();
			grad_f /= edge.length2();
			if (direction_gradient_)
			{
				grad_f -= field_grad;
			}
			grad_f *= 4 * lambda * sin(8 * angle_diff) + 2 * (1-lambda) * sin(4 * angle_diff);
		}
	}

	void DelaunayCVT::get_vertex_direction_fg(Delaunay::Vertex_handle v, double& f, vec2& grad_f) {
		f = 0.0 ;
		grad_f = vec2(0.0, 0.0);

		vec2 p0 = to_geex(v->point()) ;
		//vec2 Vd(0.0, 0.0);
		double average_length = sqrt(delaunay_->average_area_);
		double direction_factor = 0.005 * delaunay_->average_area_ * average_length * direction_factor_ ;
		Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
		Delaunay::Face_circulator jt = it ; jt++ ;
		do 
		{
			Delaunay::Vertex_handle vj = delaunay_->incident_vertices(v, jt);
			//infinite edge
			if (delaunay_->is_infinite(vj))
			{
				it++; jt++;
				continue;
			}
			//
			vec2 pj = to_geex(vj->point());
			vec2 p0pj = pj - p0;
			double E1 = 0.0;
			vec2 grad_E1(0.0, 0.0);
			double E2 = 0.0;
			vec2 grad_E2(0.0, 0.0);
			//boundary edge
			if (delaunay_->is_infinite(it) != delaunay_->is_infinite(jt))
			{
				if(direction_boundary_edge_) {
					//E1
					get_edge_direction_fg(p0pj, p0, E1, grad_E1);
					//E2...
					Delaunay::Edge_circulator e_p0pj;
					if (delaunay_->is_infinite(it))
					{
						e_p0pj = delaunay_->incident_edges(v, jt);
					}
					else {
						e_p0pj = delaunay_->incident_edges(vj, it);  //e_pjp0
					}
					CGAL::Object O_dual = delaunay_->dual(e_p0pj);
					double length_p1p2 = delaunay_->boundary_clipped_length(O_dual);
					E2 = E1 * length_p1p2;
					grad_E2 = length_p1p2 * grad_E1;
				}
			}
			//inner edge
			else {
				const vec2& p1 = it->dual ;
				const vec2& p2 = jt->dual ;
				vec2 p1p2 = p2 - p1;
				if (delaunay_->in_boundary(p1) && delaunay_->in_boundary(p2))
				{ //inner-free edge
					get_edge_direction_fg(p0pj, p0, E1, grad_E1);
					E2 = E1 * p1p2.length();
					grad_E2 = p1p2.length() * grad_E1;
				}
				else 
				{ //inner-constrained edge
					if (direction_constrained_edge_)
					{
						//E1
						get_edge_direction_fg(p0pj, p0, E1, grad_E1);
						//E2...
						Delaunay::Edge_circulator e_p0pj = delaunay_->incident_edges(v, jt);
						CGAL::Object O_dual = delaunay_->dual(e_p0pj);
						double length_p1p2 = delaunay_->boundary_clipped_length(O_dual);
						E2 = E1 * length_p1p2;
						grad_E2 = length_p1p2 * grad_E1;
					}
				}
			}

			double dir_E;
			vec2 dir_G;
			switch (direction_energy_mode_)
			{
			case ANGLE:
				dir_E = direction_factor * E1;
				dir_G = direction_factor * grad_E1;
				break;
			case ANGLE_DUAL:
				dir_E = direction_factor * E2;
				dir_G = direction_factor * grad_E2;
				break;
			case ANGLE_PRIMAL:
				dir_E = direction_factor * length(p0pj) * E1;
				dir_G = direction_factor * length(p0pj) * grad_E1;
				break;
			default:
				break;
			}

			f += dir_E;
			grad_f += dir_G;
			it++; jt++;

		} while (it != delaunay_->incident_faces(v));

		v->direction_grad = grad_f;
	}

	void DelaunayCVT::get_vertex_lloyd_fg(Delaunay::Vertex_handle v, double& f, vec2& grad_f) {
		vec2 p0 = to_geex(v->point()) ;
		vec2 Vg(0.0, 0.0) ;
		double V = 0.0 ;

		f = 0;
		if(delaunay_->dual_cell_intersects_boundary(v)) {
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			for(unsigned int i=0; i<P->size(); i++) {
				const vec2& p1 = (*P)[i].vertex[0] ;
				const vec2& p2 = (*P)[i].vertex[1] ;
				double Vi = triangle_area(p0, p1, p2) ;
				vec2 Gi = triangle_centroid(p0, p1, p2) ;
				V += Vi ;
				Vg += Vi * Gi ;
				f += Lloyd_energy(p0, p1, p2);
			}
		} else {
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			Delaunay::Face_circulator jt = it ; jt++ ;
			do {
				const vec2& p1 = it->dual ;
				const vec2& p2 = jt->dual ;
				double Vi = triangle_area(p0, p1, p2) ;
				vec2 Gi = triangle_centroid(p0, p1, p2) ;
				V += Vi ;
				Vg += Vi * Gi ;
				f += Lloyd_energy(p0, p1, p2) ;
				it++ ; jt++ ;
			} while(it != delaunay_->incident_faces(v)) ;
		}

		///dxy add
		v->area = V;
		grad_f = 2.0 * (V * p0 - Vg) ;
		v->lloyd_grad = grad_f;
	}

	bool DelaunayCVT::get_fg(Delaunay::Vertex_handle v, double& f, vec2& grad_f) {

		///dxy add
		f = 0;
		grad_f = vec2(0.0, 0.0);
		if (direction_field_)
		{
			double dir_f;
			vec2 dir_g;
			get_vertex_direction_fg(v, dir_f, dir_g);
			f += dir_f;
			grad_f += dir_g;
		}
		if (use_lloyd_energy_)
		{
			double cvt_f;
			vec2 cvt_g;
			get_vertex_lloyd_fg(v, cvt_f, cvt_g);
			f += cvt_f;
			grad_f += cvt_g;
		}
		
		bool result = true ;
		if(!delaunay_->in_boundary(to_geex(v->point()))) {
			result = false ;
		}
		return result ;
	}


	void DelaunayCVT::set_anisotropy(const vec2& X, const vec2& Y) {
		double a00 = X.x ; double a01 = Y.x ;
		double a10 = X.y ; double a11 = Y.y ;
		double d = a00 * a11 - a10 * a01 ;

		double m00 =  a11/d ; double m01 = -a01/d ;
		double m10 = -a10/d ; double m11 =  a00/d ;

		cur_func_->add_p(m00) ;
		cur_func_->add_p(m01) ;
		cur_func_->add_p(m10) ;
		cur_func_->add_p(m11) ;
	}


	void DelaunayCVT::query_anisotropy(const vec2& P, vec2& U, vec2& V) {


		switch(aniso_mode_) {
		case CONSTANT: {
			U = vec2(1.0, 0.0) ;
			V = vec2(0.0, 1.0) ;

			switch(int(X_scale_)) {
			case -2: U /= 2.0 ; break ;
			case -1: U /= ::sqrt(2.0) ; break ;
			case  0: break ;
			case  1: U *= ::sqrt(2.0) ; break ;
			case  2: U *= 2.0 ; break ;
			}

			switch(int(Y_scale_)) {
			case -2: V /= 2.0 ; break ;
			case -1: V /= ::sqrt(2.0) ; break ;
			case  0: break ;
			case  1: V *= ::sqrt(2.0) ; break ;
			case  2: V *= 2.0 ; break ;
			}
					   } break ;

		case R_INV: {
			double cx = 0.5 * (delaunay_->x_min_ + delaunay_->x_max_) ;
			double cy = 0.5 * (delaunay_->y_min_ + delaunay_->y_max_) ;

			double dx = delaunay_->x_max_ + delaunay_->x_min_ ;
			double dy = delaunay_->y_max_ + delaunay_->y_min_ ;
			double bbox_r = ::sqrt(dx * dx + dy * dy) ;
			double min_r = 0.1 * bbox_r ;

			vec2 R = P - vec2(cx, cy) ;
			double r = R.length() ;
			//            r = gx_max(r, bbox_r) ;        

			vec2 X = R ; X = normalize(X) ;
			vec2 Y(-X.y, X.x) ;


			switch(int(X_scale_)) {
			case -2: X /= r ; break ;
			case -1: X /= ::sqrt(r) ; break ;
			case  0: break ;
			case  1: X *= ::sqrt(r) ; break ;
			case  2: X *= r ; break ;
			}

			switch(int(Y_scale_)) {
			case -2: Y /= r ; break ;
			case -1: Y /= ::sqrt(r) ; break ;
			case  0: break ;
			case  1: Y *= ::sqrt(r) ; break ;
			case  2: Y *= r ; break ;
			}

			U = X ;
			V = Y ;
					} break ;

		case BORDER: {
			double theta = default_theta(P) ;
			vec2 X(cos(theta), sin(theta)) ;
			vec2 Y(-X.y, X.x) ;
			U = X ;
			V = Y ;
					 } break ;
		}

	}

	void DelaunayCVT::set_anisotropy(const vec2& P) {
		vec2 U,V ;
		query_anisotropy(P,U,V) ;
		set_anisotropy(U,V) ;
	}

	void DelaunayCVT::set_anisotropy(Delaunay::Vertex_handle V) {
		if(!use_theta_) {
			set_anisotropy(to_geex(V->point())) ;
		}
	}

	void DelaunayCVT::get_PQR(PolygonVertex& V, int center_index) {

		int nb_P = get_vertex_config(V) ;
		int nb_Q = 2 - nb_P ;

		int P[2] ;
		int Q[2] ;

		if(nb_Q == 2) {
			new_R(V) ;
		} else {
			{
				std::set<int>::iterator it = V.bisectors.begin() ;
				for(unsigned int i=0; i<nb_P; i++) {
					P[i] = *it ;
					it++ ;
				}
			}
			{
				std::set<int>::iterator it = V.boundary_edges.begin() ;
				for(unsigned int i=0; i<nb_Q; i++) {
					Q[i] = *it ;
					it++ ;
				}
			}

			new_P(delaunay_->all_vertices_[center_index]) ;
			for(unsigned int i=0; i<nb_P; i++) {
				new_P(delaunay_->all_vertices_[P[i]]) ;
			}
			for(unsigned int i=0; i<nb_Q; i++) {
				new_Q(delaunay_->boundary_[Q[i]].line()) ;
			}
		}
	}

	double DelaunayCVT::lloyd_energy() {
		double result = 0.0 ;
		for(unsigned int i=0; i<delaunay_->all_vertices_.size(); i++) {
			Delaunay::Vertex_handle v = delaunay_->all_vertices_[i] ;
			vec2 g ; double f ;
			get_fg(v, f, g) ;
			result += f ;
			v->energy = f ;
		}
		return result ;
	}

}


//------------------------- LBFGS interface ------------------------------------------------

namespace Geex {

	void DelaunayCVT::funcgrad(const double* x, double& f, double* g, bool& valid) {
		if(symbolic_) {
			funcgrad_symbolic(x, f, g, valid) ;
		} else {
			funcgrad_simple(x, f, g, valid) ;
		}
		if(has_affinity_ && affinity_weight_ > 0.0) {
			funcgrad_affinity(x,f,g) ;
		}
	}

	void DelaunayCVT::funcgrad_simple(const double* x, double& f, double* g, bool& valid) {
		valid = true ;
		set_vertices(x) ;
		std::vector<Geex::Delaunay::Vertex_handle>& all_vertices = delaunay_->all_vertices_ ;
		int cur_i = 0 ;
		f = 0.0 ;
		double gnorm2 = 0.0 ;
		for(unsigned int i=0; i<all_vertices.size(); i++) {
			double cur_f = 0.0 ; //initialization to please MSVC...
			Geex::vec2 cur_grad ;
			valid = valid && get_fg(all_vertices[i], cur_f, cur_grad) ;
			f += cur_f ;
			g[cur_i  ] = cur_grad.x ;
			g[cur_i+1] = cur_grad.y ;
			gnorm2 += cur_grad.length2() ;
			all_vertices[i]->energy = cur_f ;
			cur_i += 2 ;
			if(use_theta_) { cur_i += 2 ; }
		}
		if(lbfgs_redraw) {
			std::cerr << "Lloyd energy = " << f << std::endl ;
			std::cerr << "||g|| = " << ::sqrt(gnorm2) << std::endl ;
		}
	}

	double DelaunayCVT::default_theta(const vec2& P) {
		int edge_index = delaunay_->segments().locate(P) ;
		gx_assert(edge_index >= 0 && edge_index < delaunay_->boundary_.size()) ;
		vec2 V = delaunay_->boundary_[edge_index].vertex[1] - delaunay_->boundary_[edge_index].vertex[0] ;
		double result = 0.0 ;
		if(::fabs(V.x) < 1e-20) {
			result = M_PI / 2.0 ;
		} else {
			result = atan(V.y / V.x) ;
		}
		return result ;
	}

	double DelaunayCVT::default_theta(Delaunay::Vertex_handle v) {
		return default_theta(to_geex(v->point())) ;
		/*
		gx_assert(v->dual_intersects_boundary) ;
		Polygon2* P = delaunay_->dual_convex_clip(v, false) ;
		vec2 V(0.0, 0.0) ;
		for(unsigned int i=0; i<P->size(); i++) {
		const PolygonEdge& e = (*P)[i] ;
		V += e.vertex[1] - e.vertex[0] ;
		}
		*/
	}

	double DelaunayCVT::default_rho(Delaunay::Vertex_handle v) {

		return 1.0 ; 

		//        gx_assert(v->dual_intersects_boundary) ;        

		double cx = 0.5 * (delaunay_->x_min_ + delaunay_->x_max_) ;
		double cy = 0.5 * (delaunay_->y_min_ + delaunay_->y_max_) ;
		vec2 C(cx, cy) ;

		vec2 P = to_geex(v->point()) ;
		vec2 X = P - vec2(cx, cy) ;
		double R = X.length() ;
		return R ;
	}

	void DelaunayCVT::funcgrad_symbolic(const double* x, double& f, double* g, bool& valid) {
		valid = true ;
		set_vertices(x) ;
		std::vector<Delaunay::Vertex_handle>& all_vertices = delaunay_->all_vertices_ ;
		f = 0.0 ;
		int N = use_theta_ ? delaunay_->nb_vertices() * 4 : delaunay_->nb_vertices() * 2 ;

		for(int i=0; i<N; i++) {
			g[i] = 0.0 ;
		}
		for(unsigned int i=0; i<all_vertices.size(); i++) {
			Delaunay::Vertex_handle v = all_vertices[i] ;
			double D = 0.0 ;
			if(!delaunay_->in_boundary(to_geex(v->point()))) {
				valid = false ;
			}
			v->energy = 0.0 ;
			if(delaunay_->dual_cell_intersects_boundary(v)) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				for(unsigned int j1 = 0; j1 < P->size(); j1++) {
					int c1 = get_vertex_config((*P)[j1].vertex[0]) ;
					int c2 = get_vertex_config((*P)[j1].vertex[1]) ;

					begin_func(c1, c2) ;
					set_anisotropy(v) ;
					new_P(v) ;
					get_PQR((*P)[j1].vertex[0], v->index) ;
					get_PQR((*P)[j1].vertex[1], v->index) ;
					if(use_theta_) { 
						cur_func_->add_x(v->theta) ; 
						cur_func_->add_x(v->rho) ; 
					}
					end_func() ;
					add_to_fg(f, g) ;
					v->energy += cur_func_->f(0) ;
					D += triangle_area(to_geex(v->point()), (*P)[j1].vertex[0], (*P)[j1].vertex[1]) ;
				}
			} else {
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				Delaunay::Face_circulator jt = delaunay_->incident_faces(v) ;
				jt++ ;
				do {
					Delaunay::Vertex_handle v2 = 0 ;
					Delaunay::Vertex_handle v3 = 0 ;
					Delaunay::Vertex_handle v4 = 0 ;

					for(unsigned int j=0; j<3; j++) {
						Delaunay::Vertex_handle w1 = it->vertex(j) ;
						if(w1 != v) {
							if(jt->has_vertex(w1)) { 
								v2 = w1 ; 
							} else {
								v3 = w1 ; 
							}
						}
						Delaunay::Vertex_handle w2 = jt->vertex(j) ;
						if(w2 != v && !it->has_vertex(w2)) {
							v4 = w2 ;
						}
					}

					gx_assert(v2 != 0) ; 
					gx_assert(v3 != 0) ; 
					gx_assert(v4 != 0) ; 

					begin_func_regular() ;
					set_anisotropy(v) ;
					new_P(v) ;
					new_P(v2) ;
					new_P(v3) ;
					new_P(v4) ;
					if(use_theta_) { 
						cur_func_->add_x(v->theta) ; 
						cur_func_->add_x(v->rho) ; 
					}
					end_func() ;
					add_to_fg(f, g) ;
					v->energy += cur_func_->f(0) ;
					D += triangle_area(to_geex(v->point()), it->dual, jt->dual) ;
					it++ ;
					jt++ ;
				} while(it != delaunay_->incident_faces(v)) ;
			}
			//            v->energy /= D ;
		}
		if(lbfgs_redraw) {
			std::cerr << "Lloyd energy = " << f << std::endl ;
			double gnorm2 = 0.0 ;
			for(int i=0; i<N; i++) {
				gnorm2 += g[i]*g[i] ;
			}
			std::cerr << "||g|| = " << ::sqrt(gnorm2) << std::endl ;
		}
	}


	void DelaunayCVT::add_to_fg(double& f, double* g) {
		f += cur_func_->f(0) ;
		int nb_params = cur_P_ ;
		for(int i1 = 0; i1<nb_params; i1++) {
			if(!delaunay_->all_vertices_[iP[i1]]->locked) {
				for(int i2=0; i2<2; i2++) {
					int i = i1*2 + i2 ;
					int gi = use_theta_ ? iP[i1]*4 + i2 : iP[i1]*2 + i2 ;
					g[gi] += cur_func_->g(0,i) ;
				}
			}
		}
		if(use_theta_) {
			double dfdtheta = cur_func_->g(0,cur_func_->nb_x() - 2) ;
			double dfdrho   = cur_func_->g(0,cur_func_->nb_x() - 1) ;

			Delaunay::Vertex_handle v = delaunay_->all_vertices_[iP[0]] ;

			if(true || !v->locked) {
				if(true || !v->dual_intersects_boundary) {
					g[iP[0]*4+2] += dfdtheta ;
				}
				if(false && !v->dual_intersects_boundary) {
					g[iP[0]*4+3] += dfdrho ;
				}
			}
		}
	}

	inline bool neigh(unsigned int i, unsigned int j) {
		return (i == j+1) || (j == i+1) ;
	}

	void DelaunayCVT::load_affinity(const std::string& filename) {
		std::ifstream in(filename.c_str()) ;
		if(!in) {
			std::cerr << filename << ": could not open" << std::endl ;
			return ;
		}
		for(unsigned int i=0; i<256; i++) {
			int u1 = i % 16 ;
			int v1 = i / 16 ;
			unsigned int count = 0 ;
			for(unsigned int j=0; j<256; j++) {
				in >> affinity_[i][j] ;
				affinity_[i][j] = -affinity_[i][j]/100000.0 ;

				int u2 = j % 16 ;
				int v2 = j / 16 ;
				affinity_[i][j] = ((v1 == v2) && neigh(u1,u2) || (u1 == u2) && neigh(v1,v2)) ? -1.0 : 0.0 ;

			}
		}
		for(unsigned int i=0; i<256; i++) {
			double diag = 0 ;
			for(unsigned int j=0; j<256; j++) {
				diag += affinity_[i][j];
			}
			bool normalize = true ;
			if(normalize) {
				for(unsigned int j=0; j<256; j++) {
					affinity_[i][j] /= ::fabs(diag) ;
				}
				diag = -1.0 ;
			}
			affinity_[i][i] = -diag + 1e-20 ;
		}
		has_affinity_ = true ;
	}

	void DelaunayCVT::funcgrad_affinity(const double* x, double& f, double* g) {
		static double G[256] ;
		for(unsigned int coord = 0; coord < 2; coord++) {
			double F = 0.0 ;
			for(unsigned int i=0; i<256; i++) {
				G[i] = 0.0 ;
				for(unsigned int j=0; j<256; j++) {
					G[i] += affinity_[i][j] * x[2*j + coord] ;
				}
			}
			for(unsigned int i=0; i<256; i++) {
				F += x[2*i + coord] * G[i] ;
			} 
			std::cerr << "F_lloyd = " << f << "   F_affinity = " << F << std::endl ;

			for(unsigned int i=0; i<256; i++) {
				g[2*i+coord] += affinity_weight_ * 2.0 * G[i] ;
			}
			f += affinity_weight_ * F ;
		}
	}

}


//------------ Optimizer interface -------------

void funcgrad_cvt2d(int n, double* x, double& f, double* g) {
	bool valid ;
	Geex::DelaunayCVT* cvt = Geex::DelaunayCVT::instance() ;

	cvt->funcgrad(x,f,g, valid) ;
	if(!valid) { f += 30.0 ; }
}

void newiteration_cvt2d(int n, const double* x, double f, const double* g, double gnorm) {
	if(lbfgs_redraw) {
		glut_viewer_redraw() ;
	}
}

