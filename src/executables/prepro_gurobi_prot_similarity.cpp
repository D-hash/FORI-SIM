#define GXL_GEDLIB_SHARED

#include "src/env/ged_env.hpp"

#include <fstream>
#include <random>
#include <filesystem>
#include <chrono>
#include <cmath>
#include <regex>

#include "auxiliary/cxxopts.hpp"
#include <auxiliary/gedlib_costs.hpp>
#include <Core>
#include "auxiliary/GXLGraphReader.hpp"
#include "auxiliary/options.hpp"
#include "gurobi/FORI_VERIFICATION.hpp"
#include "auxiliary/get_lower_bound.hpp"

namespace fs = std::filesystem;


std::vector<std::string> getGXLFiles(const std::string& folderPath) {
    std::vector<std::string> gxlFiles;
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (entry.is_regular_file() && entry.path().extension() == ".gxl") {
            gxlFiles.push_back(entry.path().filename().string());
        }
    }
    return gxlFiles;
}

std::vector<std::string> selectRandomFiles(const std::vector<std::string>& files, size_t count, int seed = 0) {
    std::vector<std::string> eligibleFiles;
    for (const auto& file : files) { 
        std::string gr1 = "../data/Protein/" +  file;
        graph<std::pair<int, std::string>, std::tuple<int, int, int>> graph1 = GXLGraphReader::read_Proteins(gr1);

        eligibleFiles.push_back(file);
    }
    std::vector<std::string> shuffled = eligibleFiles;

    std::mt19937 gen(seed ? seed : std::random_device{}());

    std::shuffle(shuffled.begin(), shuffled.end(), gen);

    if (shuffled.size() > count) {
        shuffled.resize(count);
    }

    return shuffled;
}

void writeToFile(const std::vector<std::string>& files, const std::string& outputPath) {
    std::ofstream outFile(outputPath);
    if (!outFile) {
        std::cerr << "Error: Could not open file: " << outputPath << std::endl;
        return;
    }

    for (const auto& file : files) {
        outFile << file << "\n";
    }

    outFile.close();
}

int main(int argc, char **argv) {
    try {
        cxxopts::Options opts("experiment", "calc GED with given formulation, different branching priorities possible");
        opts.add_options()
            ("f, formulation", "which LP Formulation to use", cxxopts::value<std::string>())(
            "t, threads", "number of threads", cxxopts::value<int>()->default_value("1"))(
            "l, timelimit", "number of seconds the solver is allowed to run for", cxxopts::value<double>()->default_value("900"))(
            "r, seed", "random seed", cxxopts::value<int>()->default_value("1"))(
            "s, size", "number of query graphs", cxxopts::value<int>()->default_value("100"))(
            "solutionFile", "absolute path to a .gedsol file for the corresponding instance. throws runtime error if file doesn't correspond to graph1id_graph2id", cxxopts::value<std::string>()->default_value(""))(
            "v, verificationThreshold", "Set the graph similarity search threshold",cxxopts::value<double>()->default_value("1"))(
            "w, writeInFolder", "path in which to write solution file", cxxopts::value<std::string>()->default_value(""))(
            "flat, flatConstraint", "Set to 1 if model with flat constraint", cxxopts::value<bool>()->default_value("false"))(
            "p, preprocessing", "Set to 1 to use preprocessing before running gurobi", cxxopts::value<bool>()->default_value("false"))(
            "u, uniformCosts", "Set to 1 if uniform edit costs should be used, 0 otherwise", cxxopts::value<bool>()->default_value("false"));


        auto arguments = opts.parse(argc, argv);
        const std::string formulation = arguments["formulation"].as<std::string>();
        const std::string solutionFile = arguments["solutionFile"].as<std::string>();
        const std::string writePath = arguments["writeInFolder"].as<std::string>();
        const int threads = arguments["threads"].as<int>();
        const int size = arguments["size"].as<int>();
        const int seed = arguments["seed"].as<int>();
        const double timelimit = arguments["timelimit"].as<double>();
        const double threshold = arguments["verificationThreshold"].as<double>();
        const bool flat = arguments["flatConstraint"].as<bool>();
        const bool uniformCosts = arguments["uniformCosts"].as<bool>();
        const bool preproc = arguments["preprocessing"].as<bool>();
        const bool heuristic = arguments["heuristic"].as<bool>();
        options opt;
        opt.dataset_name_ = "protein";
        opt.seed_ = seed;
        opt.timelimit_ = timelimit;
        opt.threads_ = threads;
        opt.formulation_name_ = formulation;
        opt.solutionFilePath_ = solutionFile;
        opt.threshold = threshold;
        opt.flat = true;

        opt.size = size;
        opt.preprocessing_ = preproc;
        opt.heuristic_ = heuristic;
        
        if(heuristic) {
            opt.preprocessing_ = true;
        }


        std::string protfolder = "../data/Protein/";
        std::vector<std::string> protfiles = getGXLFiles(protfolder);

        std::vector<std::string> querygraphs = selectRandomFiles(protfiles, size, seed);
        
        int counter = 1;

        //setup GEDLIB env
        std::vector<int> gedlib_ids = {1, 10, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 11, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 12, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 13, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 14, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 15, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 16, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 17, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 18, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 19, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 2, 20, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 21, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 22, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 23, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 24, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 25, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 26, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 27, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 28, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 29, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 3, 30, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 31, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 32, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 33, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 34, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 35, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 36, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 37, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 38, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 39, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 4, 40, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 41, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 42, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 43, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 44, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 45, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 46, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 47, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 48, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 49, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 5, 50, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 51, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 52, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 53, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 54, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 55, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 56, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 57, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 58, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 59, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 6, 60, 600, 61, 62, 63, 64, 65, 66, 67, 68, 69, 7, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 8, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 9, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, };
        ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
        std::unordered_set<std::string> irrelevant_attributes{"aaLength"};
        std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../data/Protein", "../data/collections/allProtein.xml", ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_attributes));
        env.set_edit_costs(ged::Options::EditCosts::PROTEIN);
        env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
        env.set_method(ged::Options::GEDMethod::BRANCH);
        env.init_method();

        std::regex re(R"_(\d+)_");
        std::smatch match;


        opt.gurobi_times_.clear();
        opt.gurobi_times_.resize(querygraphs.size(), 0);
        opt.preprocessing_times_.clear();
        opt.preprocessing_times_.resize(querygraphs.size(), 0);

        for (const auto& querygraph : querygraphs) {
            std::cout << "Query graph " << counter++ << std::endl;
            std::string gr1 = "../data/Protein/" +  querygraph;
            graph<std::pair<int, std::string>, std::tuple<int, int, int>> graph1 = GXLGraphReader::read_Proteins(gr1);


            std::regex_search(querygraph, match, re);
            int id_G = std::distance(gedlib_ids.begin(), std::find(gedlib_ids.begin(), gedlib_ids.end(), std::stoi(match.str())));


            opt.Q_filename_ = graph1.get_graph_id();
            opt.Q_id_ = querygraph;
            opt.Q_num_nodes_ = std::to_string(graph1.number_of_nodes());
            opt.Q_num_edges_ = std::to_string(graph1.number_of_edges());
            opt.accepted_graphs.clear();
            opt.verification_times.clear();
            opt.upperbounds.clear();
            opt.lowerbounds.clear();
            opt.objval_ = std::numeric_limits<double>::max();
            opt.gurobi_needed.clear();
            opt.graphlist.clear();
            for (const auto& protfile : protfiles) {
                std::string gr2 = "../data/Protein-GED/Protein/" +  protfile;
                graph<std::pair<int, std::string>, std::tuple<int, int, int>> graph2 = GXLGraphReader::read_Proteins(gr2);
                opt.graphlist.push_back(protfile);


                std::regex_search(protfile, match, re);
                int id_H = std::distance(gedlib_ids.begin(), std::find(gedlib_ids.begin(), gedlib_ids.end(), std::stoi(match.str())));


                getGEDLIBcosts<std::pair<int, std::string>, std::tuple<int, int, int>> getProteinEditCosts(&graph1, &graph2);
                getProteinEditCosts.getEditCosts(uniformCosts);

                auto start = std::chrono::high_resolution_clock::now();
                if(opt.preprocessing_) {
                        env.run_method(graph_ids.at(id_G), graph_ids.at(id_H));
                        double uniform = env.get_lower_bound(graph_ids.at(id_G), graph_ids.at(id_H));

                        if(uniform > threshold) {
                            opt.objval_ = threshold + 1;
                            auto prepro_end = std::chrono::high_resolution_clock::now();
                            auto duration_prepro = std::chrono::duration_cast<std::chrono::nanoseconds>(prepro_end - start);
                            opt.preprocessing_times_[counter-2] += duration_prepro.count();
                            opt.verification_times.push_back(opt.preprocessing_times_[counter-2]);
                            continue;
                        }
                    }
                    auto end_heur = std::chrono::high_resolution_clock::now();
                    auto duration_heur = std::chrono::duration_cast<std::chrono::nanoseconds>(end_heur - start);
                    opt.preprocessing_times_[counter - 2] += duration_heur.count();

                    opt.gurobi_needed.push_back(protfile);
                    auto gurobi_start = std::chrono::high_resolution_clock::now();

                    FORI_VERIFICATION<std::pair<int, std::string>, std::tuple<int, int, int>> ilp(opt);
                    ilp.ged(graph1, graph2);

                    auto gurobi_end = std::chrono::high_resolution_clock::now();
                    auto duration_gurobi = std::chrono::duration_cast<std::chrono::nanoseconds>(gurobi_end - gurobi_start);
                    opt.gurobi_times_[counter-2] += duration_gurobi.count();
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

                    if(opt.objval_ < threshold + 1e-9) {
                        opt.accepted_graphs.push_back(protfile);
                    }
                    opt.lowerbounds.push_back(opt.final_dualbound_);
                    opt.upperbounds.push_back(opt.objval_);
                    opt.verification_times.push_back(duration.count());
                    opt.objval_ = std::numeric_limits<double>::max();
                
            }


            opt.output_ = IO::create_verification_output(opt);


            std::ostringstream stream;
            stream << std::fixed << std::setprecision(2) << threshold;
            std::string str_threshold = stream.str();
            opt.output_fname_ = writePath + "/" + opt.Q_id_ + "_" + opt.dataset_name_ + "_" + opt.formulation_name_ +
                                +"_" + std::to_string(opt.seed_) +"_" + str_threshold + (opt.flat ? "_flat" : "")
                                + (opt.preprocessing_ ? "_YESpre" : "_NOpre") + (uniformCosts ? "uniform" : "non-uniform" )+ ".json";
            IO::writeJsonToFile(opt);
        }
    }
    catch (std::exception &e) {
        std::cout << "exception " << e.what() << std::endl;
        return 1;
    }

    return 0;
}