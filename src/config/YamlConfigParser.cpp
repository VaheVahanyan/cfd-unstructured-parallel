#include "config/YamlConfigParser.hpp"

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "utils/StringUtils.hpp"

namespace
{
    template <typename T>
    void AssignIfPresent(const YAML::Node& node, const char* key, T& target) {
        if (node[key]) {
            target = node[key].as<T>();
        }
    }

    template <typename T>
    void AssignOptionalIfPresent(const YAML::Node& node, const char* key, std::optional<T>& target) {
        if (node[key]) {
            target = node[key].as<T>();
        }
    }

    [[nodiscard]] std::string ReadLowerStringRequired(const YAML::Node& node, const char* key) {
        if (!node[key]) {
            throw std::runtime_error(std::string("Missing required key: ") + key);
        }
        return utils::ToLower(node[key].as<std::string>());
    }

    [[nodiscard]] std::optional<std::string> ReadLowerStringOptional(const YAML::Node& node, const char* key) {
        if (!node[key]) {
            return std::nullopt;
        }
        return utils::ToLower(node[key].as<std::string>());
    }

    [[nodiscard]] BoundaryStateSettings ParseBoundaryState(const YAML::Node& node) {
        BoundaryStateSettings state;
        AssignIfPresent(node, "rho", state.rho);
        AssignIfPresent(node, "u", state.u);
        AssignIfPresent(node, "v", state.v);
        AssignIfPresent(node, "p", state.p);
        return state;
    }

    void SetStructuredBoundaryTag(BoundarySettings& boundary,
                                  const int tag,
                                  const std::string& type,
                                  const std::optional<BoundaryStateSettings>& state) {
        BoundaryConditionSettings bc;
        bc.type = utils::ToLower(type);
        bc.state = state;
        boundary.by_tag[tag] = bc;
    }

    [[nodiscard]] bool IsStructuredMesh(const Settings& settings) {
        return settings.mesh.source_type == MeshSourceType::StructuredCartesian;
    }

    [[nodiscard]] int EffectiveDim(const Settings& settings) {
        return settings.mesh.dim;
    }

    [[nodiscard]] StructuredMeshSettings DefaultStructuredMeshSettings() {
        return StructuredMeshSettings{};
    }
} // namespace

ParsedYamlConfig YamlConfigParser::ParseFile(const std::string& filename) const {
    const YAML::Node root = YAML::LoadFile(filename);

    ParsedYamlConfig result;

    if (root["run"]) {
        ParseRunCases(root["run"], result.run_cases);
    } else {
        result.run_cases = {"all"};
    }

    if (root["defaults"]) {
        ParseDefaults(root["defaults"], result.settings);
    }

    ValidateSettingsConsistency(result.settings);

    if (!root["cases"]) {
        throw std::runtime_error("Missing 'cases' section in YAML config");
    }

    ParseCases(root["cases"], result.initial_conditions, result.settings);
    return result;
}

void YamlConfigParser::ParseRunCases(const YAML::Node& run_node,
                                     std::vector<std::string>& run_cases) {
    run_cases.clear();

    if (!run_node["cases"]) {
        run_cases.emplace_back("all");
        return;
    }

    const YAML::Node cases_node = run_node["cases"];
    if (!cases_node.IsSequence()) {
        throw std::runtime_error("'run.cases' must be a sequence");
    }

    for (const auto& item : cases_node) {
        run_cases.push_back(item.as<std::string>());
    }
}

void YamlConfigParser::ParseDefaults(const YAML::Node& defaults_node, Settings& settings) {
    if (defaults_node["mesh"]) {
        ParseMesh(defaults_node["mesh"], settings);
    }
    if (defaults_node["physics"]) {
        ParsePhysics(defaults_node["physics"], settings);
    }
    if (defaults_node["numerics"]) {
        ParseNumerics(defaults_node["numerics"], settings);
    }
    if (defaults_node["boundary_conditions"]) {
        ParseBoundaryConditions(defaults_node["boundary_conditions"], settings);
    }
    if (defaults_node["stopping"]) {
        ParseStopping(defaults_node["stopping"], settings);
    }
    if (defaults_node["logging"]) {
        ParseLogging(defaults_node["logging"], settings);
    }
    if (defaults_node["output"]) {
        ParseOutput(defaults_node["output"], settings);
    }
    if (defaults_node["parallel"]) {
        ParseParallel(defaults_node["parallel"], settings);
    }

    ValidateSettingsConsistency(settings);
}

void YamlConfigParser::ParseCases(const YAML::Node& cases_node,
                                  std::map<std::string, InitialConditions>& initial_conditions,
                                  const Settings& defaults) {
    for (const auto& entry : cases_node) {
        const std::string case_name = entry.first.as<std::string>();
        const YAML::Node case_node = entry.second;

        InitialConditions ic;
        ic.overrides = CaseSettings{};

        ApplyCaseOverrides(case_node, ic, defaults);

        Settings effective_settings = MergeSettings(defaults, ic.overrides, CaseSettings{});
        effective_settings.simulation_case = case_name;
        ValidateSettingsConsistency(effective_settings);

        if (!case_node["initial_condition"]) {
            throw std::runtime_error("Case '" + case_name + "' does not contain 'initial_condition'");
        }

        ParseInitialCondition(case_node["initial_condition"], ic, effective_settings);
        initial_conditions[case_name] = ic;
    }
}

void YamlConfigParser::ParseMesh(const YAML::Node& node, Settings& settings) {
    const std::string mesh_type =
        node["type"] ? utils::ToLower(node["type"].as<std::string>()) : "structured";

    AssignIfPresent(node, "dim", settings.mesh.dim);

    if (mesh_type == "structured") {
        settings.mesh.source_type = MeshSourceType::StructuredCartesian;

        StructuredMeshSettings structured =
            settings.mesh.structured.has_value()
                ? *settings.mesh.structured
                : DefaultStructuredMeshSettings();

        if (node["cells"]) {
            const YAML::Node cells = node["cells"];
            AssignIfPresent(cells, "x", structured.nx);
            AssignIfPresent(cells, "y", structured.ny);
            AssignIfPresent(cells, "z", structured.nz);
        }

        if (node["domain"]) {
            const YAML::Node domain = node["domain"];

            if (domain["x_min"]) {
                structured.x_min = domain["x_min"].as<double>();
            }
            if (domain["x_max"]) {
                structured.x_max = domain["x_max"].as<double>();
            }
            if (domain["y_min"]) {
                structured.y_min = domain["y_min"].as<double>();
            }
            if (domain["y_max"]) {
                structured.y_max = domain["y_max"].as<double>();
            }

            // Backward-compatible support
            if (domain["x"]) {
                structured.x_min = 0.0;
                structured.x_max = domain["x"].as<double>();
            }
            if (domain["y"]) {
                structured.y_min = 0.0;
                structured.y_max = domain["y"].as<double>();
            }
        }

        settings.mesh.structured = structured;
        settings.mesh.gmsh_file.reset();
        return;
    }

    if (mesh_type == "gmsh_file") {
        settings.mesh.source_type = MeshSourceType::GmshFile;

        GmshFileMeshSettings gmsh_file;
        if (settings.mesh.gmsh_file) {
            gmsh_file = *settings.mesh.gmsh_file;
        }

        if (!node["source"]) {
            throw std::runtime_error("mesh.type=gmsh_file requires mesh.source");
        }

        const YAML::Node source = node["source"];
        if (!source["file"]) {
            throw std::runtime_error("mesh.source.file is required for mesh.type=gmsh_file");
        }

        gmsh_file.file_path = source["file"].as<std::string>();
        settings.mesh.gmsh_file = gmsh_file;
        settings.mesh.structured.reset();
        return;
    }

    if (mesh_type == "gmsh_geo") {
        settings.mesh.source_type = MeshSourceType::GmshGeo;

        GmshGeoMeshSettings gmsh_geo;
        if (settings.mesh.gmsh_geo) {
            gmsh_geo = *settings.mesh.gmsh_geo;
        }

        if (!node["source"]) {
            throw std::runtime_error("mesh.type=gmsh_geo requires mesh.source");
        }

        const YAML::Node source = node["source"];
        if (!source["file"]) {
            throw std::runtime_error("mesh.source.file is required for mesh.type=gmsh_geo");
        }

        gmsh_geo.file_path = source["file"].as<std::string>();
        settings.mesh.gmsh_geo = gmsh_geo;
        settings.mesh.structured.reset();
        settings.mesh.gmsh_file.reset();
        return;
    }

    if (mesh_type == "delaunay_geo") {
        settings.mesh.source_type = MeshSourceType::DelaunayGeo;

        DelaunayGeoMeshSettings delaunay_geo;
        if (settings.mesh.delaunay_geo) {
            delaunay_geo = *settings.mesh.delaunay_geo;
        }

        if (!node["source"]) {
            throw std::runtime_error("mesh.type=delaunay_geo requires mesh.source");
        }

        const YAML::Node source = node["source"];
        if (!source["file"]) {
            throw std::runtime_error("mesh.source.file is required for mesh.type=delaunay_geo");
        }

        delaunay_geo.file_path = source["file"].as<std::string>();

        settings.mesh.delaunay_geo = delaunay_geo;
        settings.mesh.structured.reset();
        settings.mesh.gmsh_file.reset();
        settings.mesh.gmsh_geo.reset();
        return;
    }

    throw std::runtime_error("Unsupported mesh.type: " + mesh_type);
}

void YamlConfigParser::ParsePhysics(const YAML::Node& node, Settings& settings) {
    if (node["eos"]) {
        settings.eos = utils::ToLower(node["eos"].as<std::string>());
    }

    if (!node["parameters"]) {
        return;
    }

    const YAML::Node parameters = node["parameters"];
    AssignIfPresent(parameters, "gamma", settings.gamma);
}

void YamlConfigParser::ParseNumerics(const YAML::Node& node, Settings& settings) {
    if (node["method"]) {
        settings.solver = utils::ToLower(node["method"].as<std::string>());
    }

    if (!node["parameters"]) {
        return;
    }

    const YAML::Node parameters = node["parameters"];

    if (parameters["time_integrator"]) {
        settings.time_integrator = utils::ToLower(parameters["time_integrator"].as<std::string>());
    }
    if (parameters["reconstruction"]) {
        settings.reconstruction = utils::ToLower(parameters["reconstruction"].as<std::string>());
    }
    if (parameters["riemann_solver"]) {
        settings.riemann_solver = utils::ToLower(parameters["riemann_solver"].as<std::string>());
    }
    if (parameters["transport_model"]) {
        settings.transport_model = utils::ToLower(parameters["transport_model"].as<std::string>());
    }

    AssignIfPresent(parameters, "cfl", settings.cfl);
    AssignIfPresent(parameters, "global_limiter", settings.global_limiter);
    AssignIfPresent(parameters, "vacuum_fix_limiter", settings.vacuum_fix_limiter);
    AssignIfPresent(parameters, "viscosity", settings.viscosity);
    AssignIfPresent(parameters, "diffusion", settings.diffusion);
}

void YamlConfigParser::ParseBoundaryConditions(const YAML::Node& node, Settings& settings) {
    BoundarySettings boundary = settings.boundary;

    if (node["tags"]) {
        const YAML::Node tags_node = node["tags"];
        if (!tags_node.IsMap()) {
            throw std::runtime_error("'boundary_conditions.tags' must be a map");
        }

        for (const auto& entry : tags_node) {
            const int tag = entry.first.as<int>();
            const YAML::Node bc_node = entry.second;

            BoundaryConditionSettings bc;
            bc.type = ReadLowerStringRequired(bc_node, "type");

            if (bc_node["state"]) {
                bc.state = ParseBoundaryState(bc_node["state"]);
            }

            boundary.by_tag[tag] = bc;
        }
    }

    // Backward-compatible structured aliases
    if (IsStructuredMesh(settings)) {
        std::optional<BoundaryStateSettings> x_min_state;
        std::optional<BoundaryStateSettings> x_max_state;
        std::optional<BoundaryStateSettings> y_min_state;
        std::optional<BoundaryStateSettings> y_max_state;
        std::optional<BoundaryStateSettings> z_min_state;
        std::optional<BoundaryStateSettings> z_max_state;

        if (node["states"]) {
            const YAML::Node states = node["states"];
            if (states["x_min"]) x_min_state = ParseBoundaryState(states["x_min"]);
            if (states["x_max"]) x_max_state = ParseBoundaryState(states["x_max"]);
            if (states["y_min"]) y_min_state = ParseBoundaryState(states["y_min"]);
            if (states["y_max"]) y_max_state = ParseBoundaryState(states["y_max"]);
            if (states["z_min"]) z_min_state = ParseBoundaryState(states["z_min"]);
            if (states["z_max"]) z_max_state = ParseBoundaryState(states["z_max"]);
        }

        if (node["default"]) {
            const std::string default_bc = utils::ToLower(node["default"].as<std::string>());

            SetStructuredBoundaryTag(boundary, k_xmin_tag, default_bc, x_min_state);
            SetStructuredBoundaryTag(boundary, k_xmax_tag, default_bc, x_max_state);

            if (settings.mesh.dim >= 2) {
                SetStructuredBoundaryTag(boundary, k_ymin_tag, default_bc, y_min_state);
                SetStructuredBoundaryTag(boundary, k_ymax_tag, default_bc, y_max_state);
            }
            if (settings.mesh.dim >= 3) {
                SetStructuredBoundaryTag(boundary, k_zmin_tag, default_bc, z_min_state);
                SetStructuredBoundaryTag(boundary, k_zmax_tag, default_bc, z_max_state);
            }
        }

        if (node["x_min"]) SetStructuredBoundaryTag(boundary, k_xmin_tag, node["x_min"].as<std::string>(), x_min_state);
        if (node["x_max"]) SetStructuredBoundaryTag(boundary, k_xmax_tag, node["x_max"].as<std::string>(), x_max_state);

        if (settings.mesh.dim >= 2) {
            if (node["y_min"])
                SetStructuredBoundaryTag(boundary, k_ymin_tag, node["y_min"].as<std::string>(),
                                         y_min_state);
            if (node["y_max"])
                SetStructuredBoundaryTag(boundary, k_ymax_tag, node["y_max"].as<std::string>(),
                                         y_max_state);
        }

        if (settings.mesh.dim >= 3) {
            if (node["z_min"])
                SetStructuredBoundaryTag(boundary, k_zmin_tag, node["z_min"].as<std::string>(),
                                         z_min_state);
            if (node["z_max"])
                SetStructuredBoundaryTag(boundary, k_zmax_tag, node["z_max"].as<std::string>(),
                                         z_max_state);
        }
    }

    settings.boundary = std::move(boundary);
}

void YamlConfigParser::ParseStopping(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "t_end", settings.t_end);
    AssignIfPresent(node, "max_steps", settings.step_end);
}

void YamlConfigParser::ParseLogging(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "every_steps", settings.log_every_steps);
    AssignIfPresent(node, "every_time", settings.log_every_time);
}

void YamlConfigParser::ParseOutput(const YAML::Node& node, Settings& settings) {
    if (node["format"]) {
        settings.output_formats = {utils::ToLower(node["format"].as<std::string>())};
    }

    if (node["formats"]) {
        if (!node["formats"].IsSequence()) {
            throw std::runtime_error("'output.formats' must be a sequence");
        }

        settings.output_formats.clear();
        for (const auto& item : node["formats"]) {
            settings.output_formats.push_back(utils::ToLower(item.as<std::string>()));
        }
    }

    AssignIfPresent(node, "directory", settings.output_dir);
    AssignIfPresent(node, "every_steps", settings.output_every_steps);
    AssignIfPresent(node, "every_time", settings.output_every_time);
}

void YamlConfigParser::ParseParallel(const YAML::Node& node, Settings& settings) {
    AssignIfPresent(node, "mpi", settings.mpi_enabled);
    AssignIfPresent(node, "omp_threads", settings.omp_threads);
    AssignIfPresent(node, "use_morton", settings.use_morton);
    if (node["decomposition"]) {
        settings.domain_decomposition_method = utils::ToLower(node["decomposition"].as<std::string>());
    }
}

void YamlConfigParser::ApplyCaseOverrides(const YAML::Node& case_node,
                                          InitialConditions& ic,
                                          const Settings& defaults) {
    CaseSettings& overrides = ic.overrides;

    if (case_node["mesh"]) {
        Settings tmp;
        ParseMesh(case_node["mesh"], tmp);
        overrides.mesh = tmp.mesh;
    }

    if (case_node["physics"]) {
        const YAML::Node physics_node = case_node["physics"];

        if (physics_node["eos"]) {
            overrides.eos = utils::ToLower(physics_node["eos"].as<std::string>());
        }

        if (physics_node["parameters"]) {
            const YAML::Node parameters = physics_node["parameters"];
            AssignOptionalIfPresent(parameters, "gamma", overrides.gamma);
        }
    }

    if (case_node["numerics"]) {
        const YAML::Node numerics_node = case_node["numerics"];

        if (numerics_node["method"]) {
            overrides.solver = utils::ToLower(numerics_node["method"].as<std::string>());
        }

        if (numerics_node["parameters"]) {
            const YAML::Node parameters = numerics_node["parameters"];

            if (parameters["time_integrator"]) {
                overrides.time_integrator = utils::ToLower(parameters["time_integrator"].as<std::string>());
            }
            if (parameters["reconstruction"]) {
                overrides.reconstruction = utils::ToLower(parameters["reconstruction"].as<std::string>());
            }
            if (parameters["riemann_solver"]) {
                overrides.riemann_solver = utils::ToLower(parameters["riemann_solver"].as<std::string>());
            }
            if (parameters["transport_model"]) {
                overrides.transport_model = utils::ToLower(parameters["transport_model"].as<std::string>());
            }

            AssignOptionalIfPresent(parameters, "cfl", overrides.cfl);
            AssignOptionalIfPresent(parameters, "global_limiter", overrides.global_limiter);
            AssignOptionalIfPresent(parameters, "vacuum_fix_limiter", overrides.vacuum_fix_limiter);
            AssignOptionalIfPresent(parameters, "viscosity", overrides.viscosity);
            AssignOptionalIfPresent(parameters, "diffusion", overrides.diffusion);
        }
    }

    if (case_node["boundary_conditions"]) {
        Settings tmp = defaults;
        if (overrides.mesh) {
            tmp.mesh = *overrides.mesh;
        }
        ParseBoundaryConditions(case_node["boundary_conditions"], tmp);
        overrides.boundary = tmp.boundary;
    }

    if (case_node["stopping"]) {
        const YAML::Node stopping_node = case_node["stopping"];
        AssignOptionalIfPresent(stopping_node, "t_end", overrides.t_end);
        AssignOptionalIfPresent(stopping_node, "max_steps", overrides.step_end);
    }

    if (case_node["logging"]) {
        const YAML::Node logging_node = case_node["logging"];
        AssignOptionalIfPresent(logging_node, "every_steps", overrides.log_every_steps);
        AssignOptionalIfPresent(logging_node, "every_time", overrides.log_every_time);
    }

    if (case_node["output"]) {
        const YAML::Node output_node = case_node["output"];

        if (output_node["format"]) {
            overrides.output_formats = std::vector<std::string>{
                utils::ToLower(output_node["format"].as<std::string>())
            };
        }

        if (output_node["formats"]) {
            if (!output_node["formats"].IsSequence()) {
                throw std::runtime_error("'output.formats' must be a sequence");
            }

            std::vector<std::string> formats;
            for (const auto& item : output_node["formats"]) {
                formats.push_back(utils::ToLower(item.as<std::string>()));
            }
            overrides.output_formats = formats;
        }

        AssignOptionalIfPresent(output_node, "directory", overrides.output_dir);
        AssignOptionalIfPresent(output_node, "every_steps", overrides.output_every_steps);
        AssignOptionalIfPresent(output_node, "every_time", overrides.output_every_time);
    }

    if (case_node["parallel"]) {
        const YAML::Node parallel_node = case_node["parallel"];
        AssignOptionalIfPresent(parallel_node, "mpi", overrides.mpi_enabled);
        AssignOptionalIfPresent(parallel_node, "use_morton", overrides.use_morton);
        AssignOptionalIfPresent(parallel_node, "omp_threads", overrides.omp_threads);
        if (parallel_node["decomposition"]) {
            overrides.domain_decomposition_method = utils::ToLower(parallel_node["decomposition"].as<std::string>());
        }
    }
}

void YamlConfigParser::ParseInitialCondition(const YAML::Node& ic_node,
                                             InitialConditions& ic,
                                             const Settings& effective_settings) {
    const std::string initial_condition_type = ReadLowerStringRequired(ic_node, "type");

    if (initial_condition_type == "structured_regions") {
        ic.type = InitialConditionType::StructuredRegions;
        ParseStructuredInitialCondition(ic_node, ic, effective_settings);
        return;
    }

    if (initial_condition_type == "constant") {
        ic.type = InitialConditionType::Constant;
        ParseConstantInitialCondition(ic_node, ic, effective_settings);
        return;
    }

    if (initial_condition_type == "region_markers") {
        ic.type = InitialConditionType::RegionMarkers;
        throw std::runtime_error("initial_condition.type=region_markers is not implemented yet");
    }

    throw std::runtime_error("Unsupported initial_condition.type: " + initial_condition_type);
}

void YamlConfigParser::ParseStructuredInitialCondition(const YAML::Node& ic_node,
                                                       InitialConditions& ic,
                                                       const Settings& effective_settings) {
    ValidateStructuredShape(ic_node, effective_settings.mesh.dim);

    StructuredRegionInitialCondition structured_ic;

    const YAML::Node interfaces = ic_node["interfaces"];
    structured_ic.interfaces_x = interfaces["x"] ? ReadVectorDouble(interfaces["x"]) : std::vector<double>{};
    structured_ic.interfaces_y = interfaces["y"] ? ReadVectorDouble(interfaces["y"]) : std::vector<double>{};

    if (effective_settings.mesh.dim == 1) {
        ParseStructured1D(ic_node, structured_ic);
    } else if (effective_settings.mesh.dim == 2) {
        ParseStructured2D(ic_node, structured_ic);
    } else {
        throw std::runtime_error("Unsupported dimension in structured initial condition");
    }

    ic.structured_regions = std::move(structured_ic);
}

void YamlConfigParser::ParseConstantInitialCondition(const YAML::Node& ic_node,
                                                     InitialConditions& ic,
                                                     const Settings& effective_settings) {
    (void)effective_settings;

    ConstantInitialCondition constant_ic;
    AssignIfPresent(ic_node, "rho", constant_ic.rho);
    AssignIfPresent(ic_node, "u", constant_ic.u);
    AssignIfPresent(ic_node, "v", constant_ic.v);
    AssignIfPresent(ic_node, "p", constant_ic.p);


    ic.constant = constant_ic;
}

void YamlConfigParser::ParseStructured1D(const YAML::Node& ic_node, StructuredRegionInitialCondition& ic) {
    const auto rho_1d = ReadVectorDouble(ic_node["rho"]);
    const auto u_1d = ReadVectorDouble(ic_node["u"]);
    const auto v_1d = ReadVectorDouble(ic_node["v"]);
    const auto p_1d = ReadVectorDouble(ic_node["p"]);

    const std::size_t nx = rho_1d.size();

    if (u_1d.size() != nx || v_1d.size() != nx || p_1d.size() != nx) {
        throw std::runtime_error("1D initial-condition arrays must have identical size");
    }

    if (nx != ic.RegionCountX() || ic.RegionCountY() != 1) {
        throw std::runtime_error("1D initial-condition shape does not match interfaces");
    }

    auto lift_1d = [](const std::vector<double>& src) -> Field2DValues {
        Field2DValues dst;
        dst.values.resize(src.size());
        for (std::size_t ix = 0; ix < src.size(); ++ix) {
            dst.values[ix].resize(1);
            dst.values[ix][0] = src[ix];
        }
        return dst;
    };

    ic.rho = lift_1d(rho_1d);
    ic.u = lift_1d(u_1d);
    ic.v = lift_1d(v_1d);
    ic.p = lift_1d(p_1d);
}

void YamlConfigParser::ParseStructured2D(const YAML::Node& ic_node, StructuredRegionInitialCondition& ic) {
    const std::size_t nx = ic.RegionCountX();
    const std::size_t ny = ic.RegionCountY();

    if (nx >= 1 && ny == 1) {
        ParseStructured1D(ic_node, ic);
        return;
    }

    if (nx == 1 && ny >= 1) {
        const auto rho_1d = ReadVectorDouble(ic_node["rho"]);
        const auto u_1d = ReadVectorDouble(ic_node["u"]);
        const auto v_1d = ReadVectorDouble(ic_node["v"]);
        const auto p_1d = ReadVectorDouble(ic_node["p"]);

        if (rho_1d.size() != ny || u_1d.size() != ny || v_1d.size() != ny ||
            p_1d.size() != ny) {
            throw std::runtime_error("2D y-only initial-condition arrays must match y-region count");
        }

        auto lift_y_only = [](const std::vector<double>& src) -> Field2DValues {
            Field2DValues dst;
            dst.values.resize(1);
            dst.values[0].resize(src.size());
            for (std::size_t iy = 0; iy < src.size(); ++iy) {
                dst.values[0][iy] = src[iy];
            }
            return dst;
        };

        ic.rho = lift_y_only(rho_1d);
        ic.u = lift_y_only(u_1d);
        ic.v = lift_y_only(v_1d);
        ic.p = lift_y_only(p_1d);

        return;
    }

    const auto rho_2d = ReadMatrixDouble(ic_node["rho"]);
    const auto u_2d = ReadMatrixDouble(ic_node["u"]);
    const auto v_2d = ReadMatrixDouble(ic_node["v"]);
    const auto p_2d = ReadMatrixDouble(ic_node["p"]);

    auto lift_2d = [&](const std::vector<std::vector<double>>& src) -> Field2DValues {
        if (src.size() != nx) {
            throw std::runtime_error("2D initial-condition x-size does not match interfaces");
        }

        Field2DValues dst;
        dst.values.resize(nx);

        for (std::size_t ix = 0; ix < nx; ++ix) {
            if (src[ix].size() != ny) {
                throw std::runtime_error("2D initial-condition y-size does not match interfaces");
            }

            dst.values[ix].resize(ny);
            for (std::size_t iy = 0; iy < ny; ++iy) {
                dst.values[ix][iy] = src[ix][iy];
            }
        }

        return dst;
    };

    ic.rho = lift_2d(rho_2d);
    ic.u = lift_2d(u_2d);
    ic.v = lift_2d(v_2d);
    ic.p = lift_2d(p_2d);
}

void YamlConfigParser::ValidateStructuredShape(const YAML::Node& ic_node, const int dim) {
    if (!ic_node["interfaces"]) {
        throw std::runtime_error("Missing 'interfaces' in structured initial condition");
    }

    if (!ic_node["rho"] || !ic_node["u"] || !ic_node["v"] || !ic_node["p"]) {
        throw std::runtime_error("structured initial condition requires rho/u/v/p");
    }

    if (dim < 1 || dim > 3) {
        throw std::runtime_error("Only dim=1,2 are supported");
    }
}

void YamlConfigParser::ValidateSettingsConsistency(const Settings& settings) {
    const int dim = settings.mesh.dim;
    if (dim < 1 || dim > 2) {
        throw std::runtime_error("mesh.dim must be 1 or 2");
    }

    if (settings.mesh.source_type == MeshSourceType::StructuredCartesian) {
        if (!settings.mesh.structured.has_value()) {
            throw std::runtime_error("Structured mesh settings are missing");
        }

        const StructuredMeshSettings& s = *settings.mesh.structured;

        if (s.nx <= 0) {
            throw std::runtime_error("structured mesh nx must be positive");
        }
        if (dim >= 2 && s.ny <= 0) {
            throw std::runtime_error("structured mesh ny must be positive for dim = 2");
        }

        if (!(s.x_max > s.x_min)) {
            throw std::runtime_error("structured mesh requires x_max > x_min");
        }
        if (dim >= 2 && !(s.y_max > s.y_min)) {
            throw std::runtime_error("structured mesh requires y_max > y_min for dim = 2");
        }
    }

    if (settings.mesh.source_type == MeshSourceType::GmshFile) {
        if (!settings.mesh.gmsh_file.has_value()) {
            throw std::runtime_error("Gmsh file mesh settings are missing");
        }
        if (settings.mesh.gmsh_file->file_path.empty()) {
            throw std::runtime_error("gmsh mesh file path is empty");
        }
    }

    if (settings.mesh.source_type == MeshSourceType::GmshGeo) {
        if (!settings.mesh.gmsh_geo.has_value()) {
            throw std::runtime_error("Gmsh geo mesh settings are missing");
        }
        if (settings.mesh.gmsh_geo->file_path.empty()) {
            throw std::runtime_error("gmsh geo file path is empty");
        }
    }

    if (settings.gamma <= 1.0) {
        throw std::runtime_error("gamma must be greater than 1");
    }
    if (settings.cfl <= 0.0) {
        throw std::runtime_error("cfl must be positive");
    }
}


std::vector<double> YamlConfigParser::ReadVectorDouble(const YAML::Node& node) {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 1D sequence");
    }

    std::vector<double> values;
    values.reserve(node.size());

    for (const auto& item : node) {
        values.push_back(item.as<double>());
    }

    return values;
}

std::vector<std::vector<double>> YamlConfigParser::ReadMatrixDouble(const YAML::Node& node) {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 2D sequence");
    }

    std::vector<std::vector<double>> values;
    values.reserve(node.size());

    for (const auto& row : node) {
        values.push_back(ReadVectorDouble(row));
    }

    return values;
}

std::vector<std::vector<std::vector<double>>> YamlConfigParser::ReadTensorDouble(const YAML::Node& node) {
    if (!node.IsSequence()) {
        throw std::runtime_error("Expected 2D sequence");
    }

    std::vector<std::vector<std::vector<double>>> values;
    values.reserve(node.size());

    for (const auto& plane : node) {
        values.push_back(ReadMatrixDouble(plane));
    }

    return values;
}

bool YamlConfigParser::HasKey(const YAML::Node& node, const char* key) {
    return static_cast<bool>(node[key]);
}
