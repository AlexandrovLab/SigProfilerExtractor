import argparse
from typing import List
from SigProfilerExtractor import sigpro


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_arguments_extractor(args: List[str], description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)

    # Core required arguments
    parser.add_argument(
        "input_type",
        help="The input file type: 'vcf', 'matrix', 'bedpe', or 'seg:TYPE'. The accepted callers for TYPE are the following {'ASCAT', 'ASCAT_NGS', 'SEQUENZA', 'ABSOLUTE', 'BATTENBERG', 'FACETS', 'PURPLE', 'TCGA'}",
    )
    parser.add_argument("output", help="Path to the output folder.")
    parser.add_argument(
        "input_data",
        help="Path to input data. For 'input_type' 'vcf' or 'bedpe', provide a path to the input folder. "
        "For 'input_type' 'matrix' or 'seg:TYPE', provide a path to the input file.",
    )

    # Optional arguments with defaults
    parser.add_argument(
        "--reference_genome",
        default="GRCh37",
        help="Reference genome (default: GRCh37).",
    )
    parser.add_argument(
        "--opportunity_genome",
        default="GRCh37",
        help="Opportunity genome (default: GRCh37).",
    )
    parser.add_argument(
        "--context_type",
        default="default",
        help="Mutational context (default: '96,DINUC,ID').",
    )
    parser.add_argument(
        "--exome",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Extract exomes (default: False).",
    )
    parser.add_argument(
        "--minimum_signatures",
        type=int,
        default=1,
        help="Minimum number of signatures (default: 1).",
    )
    parser.add_argument(
        "--maximum_signatures",
        type=int,
        default=10,
        help="Maximum number of signatures (default: 10).",
    )
    parser.add_argument(
        "--nmf_replicates", type=int, default=100, help="NMF replicates (default: 100)."
    )
    parser.add_argument(
        "--resample",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Enable resampling (default: True).",
    )
    parser.add_argument(
        "--seeds",
        default="random",
        help="Seeds for reproducible resamples, file path or 'random' (default: 'random').",
    )
    parser.add_argument(
        "--batch_size", type=int, default=1, help="Batch size for GPU (default: 1)."
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=-1,
        help="Number of processors to use (default: all available).",
    )
    parser.add_argument(
        "--gpu",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Use GPU if available (default: False).",
    )
    parser.add_argument(
        "--nmf_init",
        default="random",
        help="Initialization algorithm for NMF (default: 'random').",
    )
    parser.add_argument(
        "--precision",
        default="single",
        help="Precision for calculations (default: 'single').",
    )
    parser.add_argument(
        "--matrix_normalization",
        default="gmm",
        help="Normalization method (default: 'gmm').",
    )
    parser.add_argument(
        "--min_nmf_iterations",
        type=int,
        default=10000,
        help="Minimum NMF iterations (default: 10000).",
    )
    parser.add_argument(
        "--max_nmf_iterations",
        type=int,
        default=1000000,
        help="Maximum NMF iterations (default: 1000000).",
    )
    parser.add_argument(
        "--nmf_test_conv",
        type=int,
        default=10000,
        help="Test convergence every X iterations (default: 10000).",
    )
    parser.add_argument(
        "--nmf_tolerance",
        type=float,
        default=1e-15,
        help="NMF tolerance for convergence (default: 1e-15).",
    )
    parser.add_argument(
        "--get_all_signature_matrices",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Get all NMF matrices (default: False).",
    )
    parser.add_argument(
        "--export_probabilities",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Export probability matrix (default: True).",
    )
    parser.add_argument(
        "--stability",
        type=float,
        default=0.8,
        help="Average stability cutoff (default: 0.8).",
    )
    parser.add_argument(
        "--min_stability",
        type=float,
        default=0.2,
        help="Minimum stability cutoff (default: 0.2).",
    )
    parser.add_argument(
        "--combined_stability",
        type=float,
        default=1.0,
        help="Combined stability cutoff (default: 1.0).",
    )
    parser.add_argument(
        "--allow_stability_drop",
        type=str2bool,
        nargs="?",
        const=True,
        default=False,
        help="Allow stability drop (default: False).",
    )
    parser.add_argument(
        "--cosmic_version",
        type=float,
        default=3.4,
        help="COSMIC version for reference signatures (default: 3.4).",
    )
    parser.add_argument(
        "--make_decomposition_plots",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Generate decomposition plots (default: True).",
    )
    parser.add_argument(
        "--collapse_to_SBS96",
        type=str2bool,
        nargs="?",
        const=True,
        default=True,
        help="Collapse to SBS96 (default: True).",
    )

    return parser.parse_args(args)


class CliController:
    def dispatch_sigProfilerExtractor(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_extractor(
            user_args, "Extract mutational signatures from input samples."
        )
        sigpro.sigProfilerExtractor(
            input_type=parsed_args.input_type,
            output=parsed_args.output,
            input_data=parsed_args.input_data,
            reference_genome=parsed_args.reference_genome,
            opportunity_genome=parsed_args.opportunity_genome,
            context_type=parsed_args.context_type,
            exome=parsed_args.exome,
            minimum_signatures=parsed_args.minimum_signatures,
            maximum_signatures=parsed_args.maximum_signatures,
            nmf_replicates=parsed_args.nmf_replicates,
            resample=parsed_args.resample,
            seeds=parsed_args.seeds,
            batch_size=parsed_args.batch_size,
            cpu=parsed_args.cpu,
            gpu=parsed_args.gpu,
            nmf_init=parsed_args.nmf_init,
            precision=parsed_args.precision,
            matrix_normalization=parsed_args.matrix_normalization,
            min_nmf_iterations=parsed_args.min_nmf_iterations,
            max_nmf_iterations=parsed_args.max_nmf_iterations,
            nmf_test_conv=parsed_args.nmf_test_conv,
            nmf_tolerance=parsed_args.nmf_tolerance,
            get_all_signature_matrices=parsed_args.get_all_signature_matrices,
            export_probabilities=parsed_args.export_probabilities,
            stability=parsed_args.stability,
            min_stability=parsed_args.min_stability,
            combined_stability=parsed_args.combined_stability,
            allow_stability_drop=parsed_args.allow_stability_drop,
            cosmic_version=parsed_args.cosmic_version,
            make_decomposition_plots=parsed_args.make_decomposition_plots,
            collapse_to_SBS96=parsed_args.collapse_to_SBS96,
        )
