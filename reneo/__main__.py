"""
Entrypoint for Reneo

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    copy_config,
    echo_click,
)


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("reneo.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("reneo.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""

    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """

    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="reneo.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=8, show_default=True
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--profile",
            default=None,
            help="Snakemake profile to use",
            hidden=False,
        ),
        click.option(
            "--log",
            default="reneo.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--system_config",
            default=snake_base(os.path.join("config", "config.yaml")),
            hidden=True,
            type=click.Path(),
        ),
        click.argument("snake_args", nargs=-1),
    ]

    for option in reversed(options):
        func = option(func)
    return func


def run_options(func):
    """Reneo run-specific options"""

    options = [
        click.option(
            "--minlength",
            default=1000,
            required=False,
            help="minimum length of circular unitigs to consider",
            type=int,
            show_default=True,
        ),
        click.option(
            "--mincov",
            default=1,
            required=False,
            help="minimum coverage of paths to output",
            type=int,
            show_default=True,
        ),
        click.option(
            "--compcount",
            default=200,
            required=False,
            help="maximum unitig count to consider a component",
            type=int,
            show_default=True,
        ),
        click.option(
            "--maxpaths",
            default=10,
            required=False,
            help="maximum number of paths to resolve for a component",
            type=int,
            show_default=True,
        ),
        click.option(
            "--mgfrac",
            default=0.2,
            required=False,
            help="length threshold to consider single copy marker genes",
            type=float,
            show_default=True,
        ),
        click.option(
            "--evalue",
            default=1e-10,
            required=False,
            help="maximum e-value for vog annotations",
            type=float,
            show_default=True,
        ),
        click.option(
            "--hmmscore",
            default=50,
            required=False,
            help="minimum hmm score for vog annotations",
            type=float,
            show_default=True,
        ),
        click.option(
            "--nvogs",
            default=10,
            required=False,
            help="minimum number of vogs to consider a component",
            type=int,
            show_default=True,
        ),
        click.option(
            "--covtol",
            default=100,
            required=False,
            help="coverage tolerance for extending subpaths",
            type=int,
            show_default=True,
        ),
        click.option(
            "--alpha",
            default=1.2,
            required=False,
            help="coverage multiplier for flow interval modelling",
            type=float,
            show_default=True,
        ),
        click.option(
            "--databases",
            default=None,
            show_default=False,
            help="Custom DB location",
        ),
        click.option(
            "--hmmsearch/--no-hmmsearch",
            default=True,
            help="Perform or skip HMM searches",
            show_default=True,
        ),
        click.option(
            "--unitigs/--no-unitigs",
            default=True,
            help="Output unitigs file",
            show_default=True,
        ),
        click.option(
            "--split-paths/--no-split-paths",
            default=True,
            help="Output fasta file for each path",
            show_default=True,
        ),
    ]

    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """
    Reneo: Unraveling Viral Genomes from Metagenomes
    \b
    For more options, run:
    reneo command --help
    """
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
reneo run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           reneo run --input [assembly_graph.gfa] --reads [reads.dir/tsv]
Specify threads:    reneo run ... --threads [threads]
Disable conda:      reneo run ... --no-use-conda 
Change defaults:    reneo run ... --snake-default="-k --nolock"
Add Snakemake args: reneo run ... --dry-run --keep-going --touch
Specify targets:    reneo run ... all print_targets
Available targets:
    all             Run everything (default)
    preprocess      Run preprocessing only
    reneo           Run reneo (and preprocessing if needed)
    postprocess     Run postprocessing (with preprocessing and reneo if needed)
    print_targets   List available targets
"""


# Run command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--input",
    help="Path to assembly graph file in .GFA format",
    type=click.Path(),
    required=True,
)
@click.option(
    "--reads",
    help="Path to directory or TSV containing paired-end reads",
    type=click.Path(exists=True),
    required=True,
)
@run_options
@common_options
def run(**kwargs):
    """Run Reneo"""

    merge_config = {"reneo": kwargs}

    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "reneo.smk")),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@run_options
@common_options
def simulate(**kwargs):
    """Simulate Reneo run"""

    kwargs["input"] = snake_base(os.path.join("test_data", "assemblyGraph.gfa"))
    kwargs["reads"] = snake_base(os.path.join("test_data", "reads"))
    kwargs["databases"] = snake_base(os.path.join("test_data"))
    kwargs["snake_default"] = ["-n"]

    merge_config = {"reneo": kwargs}

    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "reneo.smk")),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@run_options
@common_options
def test(**kwargs):
    """Run Reneo with test dataset"""

    kwargs["input"] = snake_base(os.path.join("test_data", "assemblyGraph.gfa"))
    kwargs["reads"] = snake_base(os.path.join("test_data", "reads"))

    merge_config = {"reneo": kwargs}

    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "reneo.smk")),
        merge_config=merge_config,
        **kwargs
    )


# Install command
@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
def install(**kwargs):
    """Install databases"""

    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "install.smk")), **kwargs
    )


@click.command()
@common_options
def config(**kwargs):
    """Copy the system default config file"""
    copy_config(kwargs["configfile"])


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(simulate)
cli.add_command(test)
cli.add_command(install)
cli.add_command(config)
cli.add_command(citation)


def main():
    cli()


if __name__ == "__main__":
    main()
