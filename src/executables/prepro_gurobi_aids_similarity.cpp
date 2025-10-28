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
        std::string gr1 = "../data/AIDS/" +  file;
        graph<std::string, int> graph1 = GXLGraphReader::read_AIDS(gr1);

        if (graph1.number_of_nodes() >= 15) {
            eligibleFiles.push_back(file);
        }
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
        options opt;
        opt.dataset_name_ = "aids";
        opt.seed_ = seed;
        opt.timelimit_ = timelimit;
        opt.threads_ = threads;
        opt.formulation_name_ = formulation;
        opt.solutionFilePath_ = solutionFile;
        opt.threshold = threshold;
        opt.flat = true;

        opt.size = size;
        opt.preprocessing_ = preproc;
    


        std::string aidsfolder = "../data/AIDS/";
        std::vector<std::string> aidsfiles = getGXLFiles(aidsfolder);


        std::vector<std::string> querygraphs = {"20074.gxl", "42414.gxl", "33010.gxl", "27115.gxl", "435.gxl", "41217.gxl", "15750.gxl", "32612.gxl", "21643.gxl", "38188.gxl"};
        int counter = 1;

        // setup GEDLIB env
        std::vector<int> gedlib_ids = {100, 1000, 10000, 1001, 10071, 10079, 10084, 10085, 10093, 10094, 10099, 10101, 10102, 10104, 10109, 1012, 10122, 10135, 1014, 1015, 10151, 1017, 10170, 1018, 10180, 1019, 10193, 102, 10226, 10245, 10247, 10248, 10255, 10256, 10257, 1026, 10293, 10295, 1030, 10301, 10304, 10320, 10338, 1034, 1035, 1038, 10421, 1048, 10497, 10498, 1051, 10515, 10534, 10558, 10559, 10560, 10567, 10568, 10569, 1057, 10570, 10571, 10587, 10622, 10626, 10646, 1066, 10670, 10674, 10680, 107, 10706, 10725, 1074, 1079, 1080, 10830, 10839, 10841, 1085, 10858, 10868, 1088, 1089, 10890, 10894, 10908, 10928, 10934, 10936, 10946, 10960, 10962, 10963, 10966, 10971, 10972, 10973, 10979, 10982, 10992, 10999, 11016, 11017, 11018, 1102, 11020, 11036, 11041, 1105, 11059, 1106, 11064, 1108, 11092, 11095, 11099, 11100, 11108, 11130, 11146, 11147, 11166, 11168, 1118, 11182, 1119, 11200, 11215, 11216, 11227, 1127, 1132, 1136, 11369, 1137, 11377, 11385, 11386, 11390, 11392, 11394, 11396, 11423, 11424, 11425, 11426, 11447, 1150, 11509, 1152, 11521, 1154, 11563, 11564, 1157, 11573, 11575, 11576, 11582, 11586, 1159, 116, 11606, 11607, 11617, 1165, 11688, 1171, 1173, 11776, 11777, 1178, 11788, 11789, 11790, 11791, 118, 11800, 11808, 11819, 11820, 11841, 11845, 11888, 11913, 11921, 11958, 1200, 1203, 1205, 12088, 12095, 1210, 12116, 12122, 1215, 1217, 1218, 12194, 12203, 12207, 12262, 12268, 12283, 12284, 12306, 12313, 12325, 12326, 12340, 12356, 12357, 12360, 12379, 12389, 124, 12404, 12427, 12431, 12436, 12451, 12459, 1246, 12461, 12480, 12485, 12494, 125, 12519, 12535, 12537, 12538, 12541, 12542, 12544, 12555, 12556, 12560, 1259, 126, 1264, 12671, 12677, 12688, 127, 12741, 12743, 12744, 12761, 12768, 1290, 12911, 130, 1301, 13014, 13023, 13051, 13055, 13072, 1311, 1315, 1316, 1317, 1318, 1319, 13195, 132, 1320, 13210, 1322, 1323, 13231, 13232, 13233, 13239, 13240, 13270, 13272, 133, 13307, 13314, 13315, 13321, 13323, 13340, 1337, 13371, 13375, 13418, 13419, 1342, 13421, 13428, 13434, 13467, 13468, 13480, 1357, 13640, 13642, 13653, 137, 1372, 13748, 13816, 1383, 13888, 13909, 1391, 13917, 13948, 1400, 14029, 1403, 1404, 1405, 14050, 1410, 1418, 1419, 1422, 14220, 14227, 14272, 14283, 14286, 14296, 1431, 14312, 14313, 14314, 14315, 1432, 14337, 14338, 14359, 14362, 14368, 1437, 1440, 14406, 14425, 14434, 14447, 14455, 1447, 14511, 1454, 1455, 1456, 14560, 14561, 14562, 14687, 14688, 14689, 14734, 14748, 14787, 1481, 1486, 14899, 14900, 14901, 14902, 14912, 14916, 14920, 14922, 14925, 14926, 14927, 1493, 14930, 14931, 14934, 14944, 14946, 14950, 1496, 1497, 14972, 14973, 14974, 14975, 14976, 1498, 15, 150, 15065, 15069, 1511, 15135, 15210, 15224, 1525, 15283, 15286, 153, 1530, 15366, 15379, 1540, 1541, 15411, 15415, 15429, 15478, 15508, 1555, 1556, 1558, 15590, 15592, 1561, 15694, 15698, 1570, 1573, 15738, 15749, 15750, 15751, 15752, 15757, 15836, 15840, 15841, 15842, 15843, 15844, 15845, 1585, 15875, 15895, 15905, 15923, 1593, 16, 16013, 1602, 16045, 16049, 1605, 16050, 1610, 1612, 1613, 1615, 16179, 16180, 1622, 1623, 16260, 1627, 1634, 16344, 1636, 16370, 1638, 16381, 16392, 164, 1640, 16402, 1641, 165, 16582, 16586, 16592, 16598, 166, 16603, 16610, 16616, 16620, 16621, 16638, 1665, 16664, 16665, 1667, 16686, 16694, 16703, 16713, 16714, 16716, 16758, 16840, 16842, 16858, 16871, 16872, 16873, 16874, 1689, 1694, 16952, 16974, 16975, 16977, 16979, 16981, 16982, 16999, 1704, 1705, 171, 1712, 17157, 17167, 172, 1722, 17222, 17256, 1732, 1738, 1744, 17525, 17526, 17527, 17532, 17570, 17579, 17587, 17604, 17613, 17616, 17624, 17636, 1764, 1765, 17679, 1769, 1775, 1776, 17789, 17831, 17832, 17834, 17835, 1784, 17875, 17958, 17964, 1799, 18015, 18018, 18020, 18023, 18026, 18029, 18035, 18036, 18037, 18038, 18040, 18041, 1806, 18080, 18093, 18105, 1812, 1813, 18181, 1822, 1827, 1829, 18330, 18350, 18364, 18432, 18433, 18434, 18437, 18443, 18531, 1854, 1857, 1861, 18637, 18681, 18682, 187, 1877, 18788, 18791, 18797, 188, 1880, 18801, 18844, 1885, 18853, 18854, 18855, 18856, 18857, 18858, 18868, 18878, 18882, 1889, 18907, 18911, 18912, 18928, 1893, 18930, 18931, 1894, 18946, 18954, 18978, 1899, 1902, 19049, 19050, 19051, 1908, 19083, 19084, 19085, 19098, 19099, 19153, 19158, 19165, 1917, 1918, 19206, 19265, 19274, 19276, 19319, 19353, 19384, 19385, 19386, 19394, 194, 19439, 19440, 1955, 19553, 19595, 19619, 19620, 19622, 19632, 1968, 1970, 19723, 1974, 19749, 1979, 1981, 1982, 19822, 1983, 1988, 1990, 1991, 19936, 1994, 200, 2001, 20032, 2007, 20074, 20075, 20076, 20077, 20078, 2008, 20087, 20088, 20089, 2009, 20090, 201, 2010, 20115, 20116, 2013, 20140, 20200, 2025, 2028, 2029, 20366, 2040, 20416, 20417, 20418, 20425, 20449, 20453, 20478, 2051, 2052, 20556, 20566, 20643, 20646, 20680, 20681, 20682, 20683, 20684, 20685, 20688, 20691, 20692, 20693, 20696, 20697, 20698, 20700, 20702, 20706, 20708, 2075, 2076, 2078, 20782, 2080, 20805, 2097, 2098, 2099, 21020, 21044, 2105, 21072, 21080, 21086, 21087, 21147, 2117, 2118, 2119, 21251, 21274, 2128, 21290, 2130, 21332, 21333, 21371, 21376, 21377, 21378, 21379, 214, 2144, 21498, 21499, 215, 2151, 21556, 21557, 21558, 21559, 21594, 21643, 21684, 21686, 21697, 21734, 21742, 2175, 21760, 21821, 21822, 21824, 21825, 21826, 21827, 21828, 21829, 21831, 21839, 21840, 2185, 21866, 21867, 2187, 21870, 21873, 21875, 2188, 21885, 21895, 2196, 21986, 21987, 21988, 21989, 2199, 21990, 21991, 21992, 21993, 22, 22044, 2205, 22069, 22159, 22181, 22206, 22207, 22208, 22209, 22214, 2222, 2227, 22331, 2234, 22352, 22369, 2238, 2239, 2244, 22445, 22446, 22448, 2247, 2250, 22525, 22526, 22527, 22528, 22530, 22603, 22604, 22605, 22606, 2263, 22636, 22695, 22696, 22699, 22705, 22709, 22710, 22711, 22712, 22733, 2275, 22767, 2277, 22793, 228, 2291, 22994, 22995, 230, 23017, 23031, 23035, 23036, 231, 23120, 23122, 23124, 23125, 23127, 23128, 23129, 23136, 23141, 2316, 23236, 23248, 2325, 23333, 23391, 23393, 2340, 23401, 23402, 23403, 23404, 23405, 23414, 2343, 23493, 23494, 23495, 23496, 23501, 23509, 2351, 23510, 23519, 23524, 23526, 23527, 23528, 23533, 23540, 23558, 2356, 2357, 23570, 2358, 23597, 23600, 23612, 23614, 23667, 23668, 23675, 23722, 23751, 23752, 23754, 23785, 23799, 23800, 23802, 23823, 23829, 23889, 23949, 23950, 23972, 23980, 23982, 23984, 24047, 2408, 2414, 2417, 2418, 2419, 2423, 24237, 2425, 243, 244, 24418, 24419, 24420, 24421, 2445, 2449, 2454, 2458, 246, 2460, 24617, 24628, 24629, 24642, 24643, 24644, 24645, 2468, 247, 24733, 2480, 24818, 24830, 2487, 2488, 24934, 24935, 25020, 251, 252, 25256, 25279, 2528, 25316, 2532, 25353, 25354, 25355, 25357, 25358, 25359, 25360, 25367, 2548, 25529, 2553, 25530, 2554, 25558, 2563, 2566, 257, 25732, 25736, 25738, 258, 25867, 25868, 25869, 2591, 25917, 25922, 25924, 25925, 25926, 25928, 25929, 2593, 25930, 25932, 25934, 25935, 25947, 25948, 25951, 25955, 25956, 25958, 2599, 261, 26196, 262, 2623, 26266, 26284, 263, 2630, 26305, 2632, 26322, 26333, 26336, 26358, 2637, 26379, 26416, 26418, 26420, 26423, 26426, 26430, 26459, 2646, 2647, 2648, 26527, 2654, 26550, 26552, 26581, 266, 26603, 26606, 26644, 26693, 267, 26733, 26738, 26764, 26778, 26846, 2688, 269, 2693, 2694, 2696, 2697, 2709, 27115, 27116, 27117, 27136, 27137, 2715, 27162, 272, 27210, 27215, 27249, 27252, 2728, 2730, 27304, 27305, 27306, 27321, 27322, 27326, 2734, 27383, 2739, 27414, 27415, 27416, 2742, 2746, 27468, 27469, 27470, 27496, 27584, 27586, 27592, 27614, 27674, 27732, 27733, 27735, 27771, 27782, 27783, 27853, 27859, 27861, 27870, 27889, 27890, 27891, 27893, 2790, 2791, 27931, 27932, 27933, 27934, 27935, 27943, 27952, 27957, 27966, 27970, 27983, 27984, 27985, 27997, 2800, 2803, 28074, 28142, 28143, 28149, 28155, 28163, 28164, 2817, 28174, 28175, 28179, 2818, 28180, 2822, 28229, 28300, 28301, 28302, 28329, 28336, 28344, 28368, 28378, 28389, 28429, 28431, 28437, 28503, 28504, 28513, 28548, 28549, 28550, 28551, 28552, 28557, 2856, 28561, 28572, 28588, 28590, 28597, 2866, 28677, 28684, 2870, 28711, 28722, 28793, 28794, 28795, 288, 28860, 28862, 28863, 28864, 28867, 28878, 289, 28915, 2907, 29128, 29129, 29198, 292, 29385, 29387, 29433, 29437, 29440, 2947, 29473, 2948, 2949, 2950, 2951, 2952, 29523, 2954, 2955, 2956, 29586, 29598, 29646, 29648, 2965, 2966, 2967, 2968, 29690, 2970, 29731, 29736, 29782, 29783, 29785, 29790, 29812, 29816, 29845, 29886, 2995, 30, 30003, 30014, 30025, 30026, 30027, 3003, 3004, 3010, 30112, 30113, 30116, 30119, 30121, 30122, 30123, 302, 30299, 3037, 30409, 30415, 30417, 30418, 30419, 30505, 30506, 3054, 3056, 3057, 3058, 30616, 3064, 308, 30809, 30821, 3083, 30973, 31, 31026, 31079, 3108, 31080, 31081, 3109, 3110, 31110, 31111, 31112, 31113, 31114, 31137, 31148, 31256, 3128, 3141, 31438, 31439, 31440, 3147, 3149, 3151, 31510, 31511, 31569, 31570, 31574, 31631, 31634, 31639, 31640, 3180, 3181, 3182, 31897, 31898, 3190, 31997, 31998, 32167, 32228, 32304, 32363, 32488, 32611, 32612, 3266, 32692, 32701, 32791, 32862, 32867, 32872, 32873, 32875, 32903, 32904, 33010, 33062, 33064, 33065, 33085, 33103, 3311, 33111, 33112, 33115, 33132, 3315, 33169, 3317, 33170, 3324, 3325, 33403, 3344, 335, 33519, 33570, 3363, 33689, 33690, 337, 33702, 33704, 33714, 33815, 3392, 34026, 34069, 34146, 342, 34210, 34211, 34212, 34214, 34215, 34217, 34218, 34221, 34223, 34224, 34225, 34227, 34256, 34258, 34265, 34266, 34363, 3441, 34412, 34413, 34414, 34416, 3442, 34586, 34634, 347, 34730, 34759, 348, 34840, 34846, 34848, 34875, 34878, 3488, 34938, 34943, 34945, 35, 35016, 35059, 351, 35134, 35135, 35136, 35137, 35139, 35165, 35181, 35192, 35217, 35264, 35271, 3530, 35319, 35328, 35329, 35331, 35332, 35343, 35345, 35385, 35405, 35463, 35507, 35625, 357, 3577, 3578, 358, 35891, 35895, 359, 3596, 35997, 36027, 36034, 36035, 3605, 36116, 36118, 3616, 36307, 36308, 36524, 36530, 36531, 36532, 36533, 36534, 36600, 36612, 36651, 36652, 36654, 36657, 36660, 3667, 3668, 36688, 3682, 36848, 36894, 37, 3701, 37077, 37148, 37153, 3719, 37266, 37270, 37271, 37272, 373, 37312, 37328, 37508, 3754, 3759, 3760, 3761, 3762, 37633, 37634, 37660, 377, 37715, 37747, 37748, 37749, 37750, 37752, 3777, 3778, 378, 37901, 37953, 38, 380, 38056, 38188, 38257, 38311, 38312, 3833, 38436, 38585, 38632, 38639, 3883, 38872, 3892, 3894, 38940, 3895, 38951, 38965, 38966, 38970, 38971, 38974, 38975, 38982, 38985, 38986, 3905, 39056, 39262, 3937, 39372, 39374, 39375, 3938, 3940, 3943, 3945, 3951, 39511, 39527, 39570, 39692, 397, 3970, 39814, 39872, 39873, 39892, 399, 39961, 39964, 39971, 40, 400, 40102, 40117, 40201, 40261, 40265, 403, 40331, 404, 40406, 40503, 40595, 40629, 40651, 40677, 40760, 409, 40931, 4097, 41, 410, 41090, 41113, 41186, 41217, 41286, 413, 41363, 4142, 41423, 4143, 4144, 4148, 4151, 41540, 41541, 41542, 41546, 41565, 41583, 41681, 41689, 41766, 4178, 41793, 41803, 41804, 41813, 41883, 41885, 4191, 4201, 42066, 42153, 42154, 42155, 42156, 42157, 42158, 42159, 42160, 42162, 42163, 42164, 42168, 42181, 422, 42268, 42269, 42279, 42350, 42351, 42356, 42357, 42358, 42363, 42364, 42366, 42414, 4242, 42438, 4248, 4250, 4257, 426, 427, 429, 4304, 4336, 4338, 4342, 4344, 4346, 435, 4382, 4418, 4430, 4446, 4457, 4464, 4479, 449, 4493, 45, 451, 4511, 4512, 4544, 4547, 4550, 4591, 46, 461, 4614, 4619, 463, 4634, 4684, 4687, 4697, 470, 475, 4785, 4798, 4805, 4822, 4853, 486, 488, 4881, 489, 4897, 4901, 4936, 4962, 4963, 4968, 4975, 5, 5005, 5015, 5040, 5048, 5050, 5052, 5057, 5067, 5102, 52, 5209, 5253, 5257, 5268, 527, 529, 5298, 53, 5307, 5321, 5323, 5328, 5330, 5344, 537, 5372, 54, 5404, 5419, 5425, 545, 5461, 5497, 5501, 5522, 555, 5562, 5565, 5580, 5581, 5592, 5595, 56, 560, 5600, 5601, 563, 5641, 5650, 5651, 5665, 5671, 5677, 5727, 5755, 577, 5811, 5843, 59, 590, 5934, 6001, 6002, 6005, 6006, 603, 6035, 6036, 607, 6073, 6074, 6075, 6080, 6081, 6104, 6119, 6153, 6154, 6159, 6174, 6191, 62, 6204, 6217, 6223, 6227, 623, 6239, 6259, 6261, 6262, 6264, 6265, 6266, 6287, 6292, 635, 6390, 6396, 6397, 643, 6437, 6449, 645, 646, 648, 6485, 6497, 65, 66, 6614, 6616, 6643, 6656, 6678, 6687, 6691, 6701, 6733, 6749, 6757, 6786, 6787, 6801, 6851, 688, 6897, 6898, 69, 691, 6931, 695, 6958, 696, 6964, 6999, 7, 700, 7004, 7005, 7023, 7078, 7086, 7106, 7111, 7125, 718, 7191, 7204, 7217, 7280, 7281, 7319, 7320, 7340, 7400, 7411, 7414, 7438, 7440, 7444, 7469, 7470, 749, 7513, 7514, 7515, 7545, 7557, 7568, 7574, 7575, 7582, 7657, 770, 771, 7711, 7724, 7764, 7765, 7767, 7796, 7833, 7834, 7840, 7870, 7896, 7899, 79, 791, 7910, 7947, 7963, 7993, 803, 8098, 8104, 811, 8112, 8117, 8125, 813, 8130, 8153, 8162, 8165, 8201, 831, 8331, 8332, 8394, 8478, 849, 8495, 8511, 857, 8572, 8573, 8574, 8582, 859, 8608, 8628, 864, 865, 8652, 8657, 8712, 8724, 873, 8806, 881, 8890, 8894, 8914, 893, 8950, 8960, 8961, 90, 9009, 9032, 907, 91, 9101, 9103, 9105, 9130, 914, 9146, 92, 9257, 9275, 9276, 9335, 9345, 9346, 9367, 9368, 9386, 94, 941, 948, 95, 950, 9600, 961, 9615, 9642, 9680, 9690, 9693, 97, 973, 975, 9770, 9796, 9797, 983, 9834, 986, 987, 9901, 9902, 9903, 991, 9915, 9927, 9931, 9939, 9945, 9949, 9951, 9982, 9987};
        ged::GEDEnv<ged::GXLNodeID, ged::GXLLabel, ged::GXLLabel> env;
        std::unordered_set<std::string> irrelevant_attributes{"x", "y", "symbol", "charge"};
        std::vector<ged::GEDGraph::GraphID> graph_ids(env.load_gxl_graphs("../data/AIDS/", "../data/collections/allAids.xml", ged::Options::GXLNodeEdgeType::LABELED, ged::Options::GXLNodeEdgeType::LABELED, irrelevant_attributes));
        if(uniformCosts){
            env.set_edit_costs(ged::Options::EditCosts::CONSTANT);
            env.set_method(ged::Options::GEDMethod::BRANCH_UNIFORM);
        }
        else{
            env.set_edit_costs(ged::Options::EditCosts::CHEM_2);
            env.set_method(ged::Options::GEDMethod::BRANCH);
        }
        env.init(ged::Options::InitType::EAGER_WITHOUT_SHUFFLED_COPIES);
        env.init_method();

        std::regex re(R"_(\d+)_");
        std::smatch match;


        opt.gurobi_times_.clear();
        opt.gurobi_times_.resize(querygraphs.size(), 0);
        opt.preprocessing_times_.clear();
        opt.preprocessing_times_.resize(querygraphs.size(), 0);


        for (const auto& querygraph : querygraphs) {
            std::cout << "Query graph " << counter++ << std::endl;
            std::string gr1 = "../data/AIDS/" +  querygraph;
            graph<std::string, int> graph1 = GXLGraphReader::read_AIDS(gr1);


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
            for (const auto& aidsfile : aidsfiles) {
                std::string gr2 = "../data/AIDS/" +  aidsfile;
                graph<std::string, int> graph2 = GXLGraphReader::read_AIDS(gr2);
                opt.graphlist.push_back(aidsfile);


                std::regex_search(aidsfile, match, re);
                int id_H = std::distance(gedlib_ids.begin(), std::find(gedlib_ids.begin(), gedlib_ids.end(), std::stoi(match.str())));


                getGEDLIBcosts<std::string, int> getAIDSEditCosts(&graph1, &graph2);
                getAIDSEditCosts.getEditCosts(uniformCosts);

                auto start = std::chrono::high_resolution_clock::now();
                unsigned int lb = compute_lower_bound(graph1, graph2, threshold);
                if(lb <= threshold) {
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

                    opt.gurobi_needed.push_back(aidsfile);
                    auto gurobi_start = std::chrono::high_resolution_clock::now();

                    FORI_VERIFICATION<std::string, int> ilp(opt);
                    ilp.ged(graph1, graph2);

                    auto gurobi_end = std::chrono::high_resolution_clock::now();
                    auto duration_gurobi = std::chrono::duration_cast<std::chrono::nanoseconds>(gurobi_end - gurobi_start);
                    opt.gurobi_times_[counter-2] += duration_gurobi.count();
                    auto end = std::chrono::high_resolution_clock::now();
                    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

                    if(opt.objval_ < threshold + 1e-9) {
                        opt.accepted_graphs.push_back(aidsfile);
                    }
                    opt.lowerbounds.push_back(opt.final_dualbound_);
                    opt.upperbounds.push_back(opt.objval_);
                    opt.verification_times.push_back(duration.count());
                    opt.objval_ = std::numeric_limits<double>::max();
                } else {
                    auto lb_compute_end = std::chrono::high_resolution_clock::now();
                    auto duration_lb_compute = std::chrono::duration_cast<std::chrono::nanoseconds>(lb_compute_end - start);
                    opt.preprocessing_times_[counter - 2] = duration_lb_compute.count();
                    opt.verification_times.push_back(opt.preprocessing_times_[counter-2]);
                }
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