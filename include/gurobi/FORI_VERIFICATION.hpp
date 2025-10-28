#ifndef GEDC_FORI_VERIFICATION_HPP
#define GEDC_FORI_VERIFICATION_HPP

#include "gurobi_c++.h"
#include "auxiliary/graph.hpp"
#include "utils.hpp"

#include "auxiliary/io.hpp"
#include "auxiliary/options.hpp"

template<typename T, typename U>
class FORI_VERIFICATION {
private:

        bool testing_ = false;
        double constant_ = 0.0;
        bool relax_ = false;
        char VAR_TYPE_ = GRB_BINARY;
        int threads_ = 1;
        int seed_ = 1;

public:

        FORI_VERIFICATION() = default;

        explicit FORI_VERIFICATION(options &opt) : opt_(opt) {}


        inline void relax() {
                relax_ = true;
                opt_.relaxed_ = true;
        }

        inline void set_testing() { testing_ = true; }

        inline void set_threads(int threads) { threads_ = threads; }

        inline void set_seed(int seed) { seed_ = seed; }

        std::string log_name_;
        std::string output_fname_;
        double timelimit_ = 0;
        options &opt_;



        std::string ged(graph<T, U> &G, graph<T, U> &H) {

                auto n_g = G.number_of_nodes(); 
                auto n_h = H.number_of_nodes(); 
                auto m_g = G.number_of_edges(); 
                auto m_h = H.number_of_edges();

                try {
                        GRBEnv *env;
                        GRBVar **node_sub = new GRBVar*[n_g];
                        for (int i = 0; i < n_g; i++) {
                                node_sub[i] = new GRBVar[n_h];
                        }
                        GRBVar **edge_sub = new GRBVar*[m_g];
                        for (int i = 0; i < m_g; i++) {
                                edge_sub[i] = new GRBVar[m_h];
                        }
                        GRBVar **edge_sub_rev = new GRBVar*[m_g];
                        for (int i = 0; i < m_g; i++) {
                                edge_sub_rev[i] = new GRBVar[m_h];
                        }

                        env = new GRBEnv();
                        GRBModel *model = new GRBModel(*env);
                        model->set(GRB_StringAttr_ModelName, "FORI_VERIFICATION_"+opt_.dataset_name_ + "_"+opt_.G_id_ + "_" +opt_.H_id_);
                        model->set(GRB_StringParam_LogFile, opt_.log_fname_ + ".log");
                        model->set(GRB_IntParam_Threads, opt_.threads_);

                        model->set(GRB_DoubleParam_TimeLimit, opt_.timelimit_);


                        model->set(GRB_IntParam_Seed, opt_.seed_);

                        model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

                        std::vector<std::vector<double>> c_ik(n_g, std::vector<double>(n_h, 0));
                        std::vector<double> c_ie(n_g, 0);
                        std::vector<double> c_ek(n_h, 0);
                        std::vector<std::vector<double>> c_ijkl(m_g, std::vector<double>(m_h, 0));
                        std::vector<double> c_ije(m_g, 0);
                        std::vector<double> c_ekl(m_h, 0);
                        GRBLinExpr objfunc;

                        G.cost_function(G, H, c_ik, c_ie, c_ek, c_ijkl, c_ije, c_ekl);

                        if (G.get_dataset() == "CMU") {
                                for (int i = 0; i < n_g; i++) {
                                        c_ie[i] = 0;
                                }
                                for (int k = 0; k < n_h; k++) {
                                        c_ek[k] = 0;
                                }
                        }

                        for (int i = 0; i < n_g; i++) {
                                for (int k = 0; k < n_h; k++) {
                                        node_sub[i][k] = model->addVar(0, 1, 0, GRB_BINARY,
                                                                       "x" + std::to_string(i) + "_" + std::to_string(k));
                                        objfunc += (c_ik[i][k] - c_ie[i] - c_ek[k]) * node_sub[i][k];
                                }
                        }


                        for (int ij = 0; ij < m_g; ij++) {
                                for (int kl = 0; kl < m_h; kl++) {
                                        edge_sub[ij][kl] = model->addVar(0, 1, 0,
                                                                         GRB_BINARY, "z_" + std::to_string(G.get_edge(ij).first) + "_" +
                                                                                         std::to_string(G.get_edge(ij).second) + "_" + std::to_string(H.get_edge(kl).first) + "_" +
                                                                                         std::to_string(H.get_edge(kl).second));
                                        objfunc += (c_ijkl[ij][kl] - c_ije[ij] - c_ekl[kl]) * edge_sub[ij][kl];
                                        edge_sub_rev[ij][kl] = model->addVar(0, 1, 0,
                                                                 GRB_BINARY, "z_" +  std::to_string(G.get_edge(ij).first) + "_" +
                                                                                 std::to_string(G.get_edge(ij).second) + "_" + std::to_string(H.get_edge(kl).second) + "_" +
                                                                                 std::to_string(H.get_edge(kl).first));
                                        objfunc += (c_ijkl[ij][kl] - c_ije[ij] - c_ekl[kl]) * edge_sub_rev[ij][kl];
                                }
                        }


                        GRBLinExpr le;
                        for (int i = 0; i < n_g; i++) {
                                le = 0;
                                for (int k = 0; k < n_h; k++)
                                        le += node_sub[i][k];
                                model->addConstr(le, GRB_LESS_EQUAL, 1, "Ass_G_"+std::to_string(i));
                        }

                        GRBLinExpr le1;
                        for (int k = 0; k < n_h; k++) {
                                le1 = 0;
                                for (int i = 0; i < n_g; i++)
                                        le1 += node_sub[i][k];
                                model->addConstr(le1, GRB_LESS_EQUAL, 1, "Ass_H_"+std::to_string(k));
                        }

                        GRBLinExpr le2;
                        GRBLinExpr le3;
                        for (int ij = 0; ij < m_g; ij++) {
                                for (int k = 0; k < n_h; k++) {
                                        le2 = 0;
                                        le3 = node_sub[G.get_edge(ij).first][k];
                                        for (int kl = 0; kl < m_h; kl++) {
                                                if (k == H.get_edge(kl).first) {
                                                        le2 += edge_sub[ij][kl];
                                                }
                                                else if(k == H.get_edge(kl).second) {
                                                          le2 += edge_sub_rev[ij][kl];
                                                }
                                        }
                                        model->addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_1_G(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                   std::to_string(G.get_edge(ij).second) + ")" + "_" + std::to_string(k));
                                }
                        }

                        le2.clear();
                        le3.clear();
                        for (int ij = 0; ij < m_g; ij++) {
                                for (int k = 0; k < n_h; k++) {
                                        le2.clear();
                                        le3.clear();
                                        le3 = node_sub[G.get_edge(ij).second][k];
                                        for (int kl = 0; kl < m_h; kl++) {
                                                if (k == H.get_edge(kl).first) {
                                                        le2 += edge_sub_rev[ij][kl];
                                                }
                                                else if(k == H.get_edge(kl).second) {
                                                        le2 += edge_sub[ij][kl];
                                                }
                                        }
                                        model->addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_2_G_(" + std::to_string(G.get_edge(ij).first) + "," +
                                                                                   std::to_string(G.get_edge(ij).second) + ")" + "_" + std::to_string(k));
                                }
                        }

                        le2.clear();
                        le3.clear();

                        for (int kl = 0; kl < m_h; kl++) {
                          for (int i = 0; i < n_g; i++) {
                                  le2.clear();
                                  le3.clear();
                            le3 = node_sub[i][H.get_edge(kl).first];
                            for (int ij = 0; ij < m_g; ij++) {
                              if (i == G.get_edge(ij).first) {
                                le2 += edge_sub[ij][kl];
                              }
                              if (i == G.get_edge(ij).second) {
                                le2 += edge_sub_rev[ij][kl];
                              }
                            }
                                  model->addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_H_(" + std::to_string(H.get_edge(kl).first) + "," +
                                            std::to_string(H.get_edge(kl).second) + ")" + "_" +
                                            std::to_string(i));
                          }
                        }

                        le2.clear();
                        le3.clear();

                        for (int kl = 0; kl < m_h; kl++) {
                                for (int i = 0; i < n_g; i++) {
                                        le2.clear();
                                        le3.clear();
                                        le3 = node_sub[i][H.get_edge(kl).second];
                                        for (int ij = 0; ij < m_g; ij++) {
                                                if (i == G.get_edge(ij).first) {
                                                        le2 += edge_sub_rev[ij][kl];
                                                }
                                                if (i == G.get_edge(ij).second) {
                                                        le2 += edge_sub[ij][kl];
                                                }
                                        }
                                        model->addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_H_(" + std::to_string(H.get_edge(kl).second) + "," +
                                                  std::to_string(H.get_edge(kl).first) + ")" + "_" +
                                                  std::to_string(i));
                                }
                        }

                        for (int i = 0; i < n_g; i++)
                                constant_ += c_ie[i];

                        for (int k = 0; k < n_h; k++)
                                constant_ += c_ek[k];

                        for (int ij = 0; ij < m_g; ij++)
                                constant_ += c_ije[ij];

                        for (int kl = 0; kl < m_h; kl++)
                                constant_ += c_ekl[kl];

                        auto cons = model->addVar(1.0,1.0,0,GRB_BINARY,"constant");
                        objfunc += cons*constant_;
                        model->update();
                        model->setObjective(objfunc);
                        model->set(GRB_IntParam_OutputFlag, 1);

                        for (int i = 0; i < model->get(GRB_IntAttr_NumVars); ++i) {
                                GRBVar var = model->getVar(i);
                                var.set(GRB_CharAttr_VType, GRB_CONTINUOUS);
                        }
                        model->update();
                        model->set(GRB_IntParam_Method, 2);
                        
                        model->set(GRB_IntParam_Crossover, 0);
                        model->optimize();
                        if(model->get(GRB_DoubleAttr_ObjVal) < opt_.threshold + 1e-6){
                                model->reset();
                                model->addConstr(objfunc <= opt_.threshold, "threshold");
                                for (int i = 0; i < model->get(GRB_IntAttr_NumVars); ++i) {
                                        GRBVar var = model->getVar(i);
                                        var.set(GRB_CharAttr_VType, GRB_BINARY);
                                }
                                model->update();
                                model->set(GRB_IntParam_Method, 2);
                        
                                model->set(GRB_IntParam_Crossover, 0);
                                model->set(GRB_IntParam_SolutionLimit, 1);
                        model->addConstr(objfunc <= opt_.threshold, "threshold");
                        model->set(GRB_IntParam_SolutionLimit, 1);
                        model->optimize();
                        }
                        else {
                                opt_.objval_ = model->get(GRB_DoubleAttr_ObjVal);
                                opt_.final_dualbound_ = model->get(GRB_DoubleAttr_ObjVal);
                                delete[] node_sub;
                                delete[] edge_sub;
                                delete[] edge_sub_rev;
                                delete model;
                                delete env;

                                return "OK";
                        }
                                
                        
                        model->optimize();

                      
                        opt_.constant_ = constant_;

                        if(opt_.flat) {
                                int status = model->get(GRB_IntAttr_Status);

                                if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE) {
                                        opt_.objval_ = std::numeric_limits<double>::max();
                                        opt_.final_dualbound_ = 0;
                                        delete[] node_sub;
                                        delete[] edge_sub;
                                        delete[] edge_sub_rev;
                                        delete model;
                                        delete env;

                                        return "OK";
                                }
                        }
                        opt_.bbnodecount_ = model->get(GRB_DoubleAttr_NodeCount);
                        opt_.mipgap_ = model->get(GRB_DoubleAttr_MIPGap);
                        opt_.objval_ = model->get(GRB_DoubleAttr_ObjVal);
                        opt_.final_dualbound_ = model->get(GRB_DoubleAttr_ObjBound);
                        opt_.time_ = model->get(GRB_DoubleAttr_Runtime);
                        opt_.status_ = model->get(GRB_IntAttr_Status);
                        opt_.n_vars_full_ = model->get(GRB_IntAttr_NumVars);
                        opt_.n_cons_full_ = model->get(GRB_IntAttr_NumConstrs);
                        opt_.numNZ_ = model->get(GRB_IntAttr_NumNZs);
                        opt_.iterCount_ = model->get(GRB_DoubleAttr_IterCount);

                        delete[] node_sub;
                        delete[] edge_sub;
                        delete[] edge_sub_rev;
                        delete model;
                        delete env;

                        return "OK";
                }
                catch (GRBException &e) {
                        std::cout << "Error code = " << e.getErrorCode() << std::endl;
                        std::cout << e.getMessage() << std::endl;
                }
                exit(1);
        }

};


#endif //GEDC_FORI_VERIFICATION_HPP