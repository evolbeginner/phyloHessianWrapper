#! /bin/env ruby


#########################################
# last updated 2025-08-21
# well logged


#########################################
$VERBOSE = nil
DIR = File.dirname($0)
ADD_SCRIPTS_DIR = File.join(DIR, "additional_scripts")
$VERBOSE = true


#########################################
require 'getoptlong'
require 'tempfile'
require 'fileutils'
require 'colorize'

require_relative 'additional_scripts/Dir'
require_relative 'additional_scripts/util'


#########################################
JULIA_BL = File.join(DIR, "julia_bl.jl")
GEN_BASICS = File.join(DIR, 'parse_tree_seq.R')
GEN_BRANCH = File.join(ADD_SCRIPTS_DIR, 'generate_branch_out_mat.sh')

IQTREE = 'iqtree'
PHYML = 'phyml'
JULIA = 'julia'
RUBY = 'ruby'

REGULAR_DIR = File.join(DIR, 'substitution_model', 'regular')
MFM_DIR = File.join(DIR, 'substitution_model', 'mfm')

DO_MCMCTREE = File.expand_path("~/lab-tools/dating/do_mcmctree.rb")
BS_INBV = File.expand_path("~/lab-tools/dating/hessian/create_hessian_by_bootstrapping.rb")
MFATOPHY = File.expand_path(File.join(ADD_SCRIPTS_DIR, 'MFAtoPHY.jl'))


#########################################
# Initialize logging
def init_logging(outdir)
  timestamp = Time.now.strftime("%Y%m%d_%H%M%S")
  log_file = File.join(outdir, "execution_#{timestamp}.log")
  $log = File.open(log_file, 'w')
  $log.sync = true # Make sure logs are written immediately
end

def log(message)
  #puts message
  $log.puts message
end

def log_error(message)
  STDERR.puts message.red
  $log.puts "ERROR: #{message}"
end

def log_command(cmd)
  log("Executing: #{cmd}")
end

def execute_command(cmd, cmd_stdout=nil)
  log_command(cmd)
  STDOUT.puts cmd_stdout unless cmd_stdout.nil?
  output = `#{cmd} 2>&1`
  if $?.success?
    log(output) unless output.empty?
  else
    log_error("Command failed with status #{$?.exitstatus}: #{cmd}")
    log_error(output)
    raise "Command execution failed"
  end
  output
end


#########################################
def check_model_num(a, model)
  if a.size == 1
    return a[0]
  elsif a.size >= 2
    raise "Error! More than one mfm is specified in #{model}"
  else
    return 'nothing'
  end
end


def parse_model_name(model)
  non_mfms = []
  mfms = []

  model.split('+').each do |m|
    if %w[dat nex].any? { |suffix| File.exist?(File.join(REGULAR_DIR, m + '.' + suffix)) }
      non_mfms << m
    elsif %w[dat nex].any? { |suffix| File.exist?(File.join(MFM_DIR, m + '.' + suffix)) }
      mfms << m
    end
  end

  non_mfm = check_model_num(non_mfms, model)
  non_mfm = non_mfm == 'nothing' ? 'POISSON' : non_mfm
  mfm = check_model_num(mfms, model)

  [non_mfm, mfm]
end


def run_mcmctree(outdirs, tree_indir, phylip, clock, bd, bsn)
  outdirs[:date] = File.join(File.dirname(outdirs[:inBV]), 'date')
  mkdir_with_force(outdirs[:date])
  inBV_file = File.join(outdirs[:inBV], 'in.BV')
  add_argu = "--inBV #{inBV_file}"
  cmd = "#{RUBY} #{DO_MCMCTREE} --outdir #{outdirs[:date]} --tree_indir #{tree_indir} -i #{phylip} --prot --clock #{clock} --bd #{bd} --bsn #{bsn} --print 2 #{add_argu} --force --one_codeml"
  execute_command(cmd, "running MCMCtree ......")
end


def create_phylip(seqfile)
  tmp_dir = Dir.mktmpdir
  cmd = "julia #{MFATOPHY} #{seqfile}"
  execute_command(cmd)
  FileUtils.mv("#{seqfile}.phy", tmp_dir)
  File.join(tmp_dir, File.basename(seqfile) + '.phy')
end


def run_iqtree(treefile, model, pmsf, seqfile, outdir, blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
  output_file = File.join(outdir, 'iqtree_output.log')
  
  if !pmsf.nil?
    cmd = "#{IQTREE} -te #{treefile} -m #{model} -s #{seqfile} -pre #{outdir}/iqtree -redo -wsl -quiet #{iqtree_add_arg0} -T #{cpu}"
  end
  
  cmd = "#{IQTREE} -blfix #{treefile} -m POISSON -s #{seqfile} -pre #{outdir}/iqtree -redo -wsl -quiet #{iqtree_add_arg0} -T #{cpu}"
  ` #{cmd} `
  
  cmd = "#{IQTREE} -te #{outdir}/iqtree.treefile -m #{model} -s #{seqfile} -pre #{outdir}/iqtree -redo -wsl -quiet #{iqtree_add_arg} -blmin #{blmin} -T #{cpu}"
  execute_command(cmd, "Running IQ-Tree ......")
end


def run_phyml(treefile, model, pmsf, seqfile, outdir, blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
  if pmsf.nil?
    log_error("-r ref_treefile not given! Exiting ......")
    exit(1)
  end
  
  phylip = create_phylip(seqfile)
  FileUtils.cp(phylip, File.join('.', "#{getCorename(seqfile)}.phy"))
  phylip = "#{getCorename(seqfile)}.phy"

  phyml_model_arg = iqtree_to_phyml(model)
  cmd = "#{PHYML} -i #{phylip} -u #{treefile} -d aa -o lr #{phyml_model_arg}"
  
  log_file = File.join(outdir, 'phyml_output.log')
  begin
    execute_command("#{cmd} > #{log_file} 2>&1", "Running PhyML ......")
  rescue
    cmd2 = "phyml-mixSS -i #{phylip} -u #{treefile} -d aa -o lr #{phyml_model_arg}"
    log("First PhyML command failed, trying alternative: #{cmd2}")
    execute_command("#{cmd2} > #{log_file} 2>&1", "Running phyml-mixSS ......")
  end
end


#########################################
def show_help
  puts <<~HELP
    #{File.basename($0)} - v0.2

    Usage: #{File.basename($0)} [options]

    Required Options:
      -s, --aln FILE            Input sequence alignment file
      -t, --tree FILE           Input tree file
      -r, --ref_tree FILE       Reference tree file
      --outdir DIR              Output directory

    Analysis Options:
      -m, --model MODEL         Substitution model (default: POISSON)
      --mfm MODEL               Mixture model specification
      --st TYPE                 Sequence type (AA/NT, default: AA)
      --sf FORMAT               Sequence format (fasta/phylip, default: fasta)
      --pmsf TYPE               PMSF type (relaxed/stringent)
      --phylo_prog PROGRAM      Phylogeny program (iqtree/phyml, default: iqtree)
      --cpu N                   Number of CPUs (default: 4)

    Branch Length Options:
      --blmin FLOAT             Minimum branch length (default: 4e-6)
      --hessian_type TYPE       Hessian calculation method (default: SKT2004)

    MCMCTree Options:
      --run_mcmctree            Run MCMCTree after analysis
      --clock MODEL             Clock model (default: IR)
      --bd RATES                Birth-death rates (default: 1,1,0.1)
      --bsn BRANCHES            Branches for sigma2 (default: 1,2,3)

    Other Options:
      --force                   Overwrite existing files
      -h, --help                Show this help message

    Examples:
      only in.BV:
        #{File.basename($0)} -s alignment.fa -t tree.nwk -m LG+C60+G -r ref_tree.nwk --outdir only_inBV

      Run full pipeline:
        #{File.basename($0)} -s alignment.fa -t tree.nwk -r ref_tree.nwk --run_mcmctree --outdir full_analysis

      Citation:
        Wang S, Meade A. Molecular Clock Dating of Deep-Time Evolution Using Complex Mixture Models. https://www.biorxiv.org/content/10.1101/2025.07.17.665246v1. bioRxiv. 2025:2025-07.
  HELP
  exit
end


#########################################
#########################################
seqfile = nil
st = 'AA'
seq_format = 'fasta'
treefile = nil
model = 'POISSON'
mfm = 'nothing'
ref_treefile = nil
phylo_prog = :iqtree
pmsf = nil
is_mwopt = true
cpu = 4
outdir = nil
is_force = false

is_run_mcmctree = false
tree_indir = nil
phylip = nil
clock = 'IR'
bd = '1,1,0'
bsn = '1,2,3'

blmin = 4e-6
iqtree_add_arg0 = ['-mwopt', '-keep-ident'].join(' ')
iqtree_add_arg = iqtree_add_arg0
julia_bl_add_arg = nil

hessian_type = 'STK2004'
iqtree_indir = nil

outdirs = Hash.new


#########################################
opts = GetoptLong.new(
  ['-s', '--aln', GetoptLong::REQUIRED_ARGUMENT],
  ['--st', GetoptLong::REQUIRED_ARGUMENT],
  ['--sf', '-f', GetoptLong::REQUIRED_ARGUMENT],
  ['-t', '--tree', GetoptLong::REQUIRED_ARGUMENT],
  ['-m', '--model', GetoptLong::REQUIRED_ARGUMENT],
  ['--mfm', GetoptLong::REQUIRED_ARGUMENT],
  ['-r', '--reftre', '--ref_tre', '--reftree', '--ref_tree', GetoptLong::REQUIRED_ARGUMENT],
  ['--cpu', GetoptLong::REQUIRED_ARGUMENT],
  ['--pmsf', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--run_mcmctree', GetoptLong::NO_ARGUMENT],
  ['--tree_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--phylip', GetoptLong::REQUIRED_ARGUMENT],
  ['--clock', GetoptLong::REQUIRED_ARGUMENT],
  ['--bd', GetoptLong::REQUIRED_ARGUMENT],
  ['--bsn', GetoptLong::REQUIRED_ARGUMENT],
  ['--hessian_type', GetoptLong::REQUIRED_ARGUMENT],
  ['--no_mwopt', GetoptLong::NO_ARGUMENT],
  ['--iqtree_indir', GetoptLong::REQUIRED_ARGUMENT],
  ['--phylo_prog', GetoptLong::REQUIRED_ARGUMENT],
  ['--phyml', GetoptLong::NO_ARGUMENT],
  ['--iqtree', GetoptLong::NO_ARGUMENT],
  ['-h', GetoptLong::NO_ARGUMENT]
)

opts.each do |opt, value|
  case opt
    when '-s', '--aln'
      seqfile = value
    when '--st'
      seq_type = value
    when '-f', '--sf'
      seq_format = value
    when '-t'
      treefile = value
    when '-m', '--model'
      model = value
    when '--mfm'
      mfm = value
    when '--cpu'
      cpu = value.to_i
    when '--pmsf'
      pmsf = value
    when '-r', '--reftre', '--ref_tre', '--ref_tree', '--reftree'
      ref_treefile = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--run_mcmctree'
      is_run_mcmctree = true
    when '--tree_indir'
      tree_indir = value
    when '--phylip'
      phylip = value
    when '--clock'
      clock = value
    when '--bd'
      bd = value
    when '--bsn'
      bsn = value
    when '--hessian_type'
      hessian_type = value
    when '--no_mwopt'
      is_mwopt = false
    when '--phylo_prog'
      phylo_prog = value.to_sym
    when '--phyml'
      phylo_prog = :phyml
    when '--iqtree'
      phylo_prog = :iqtree
    when '--iqtree_indir'
      iqtree_indir = value
    when '-h'
      show_help()
  end
end

if not is_mwopt
  iqtree_add_arg0.gsub!(' -mwopt', '')
  iqtree_add_arg.gsub!(' -mwopt', '')
end

model = model.dup
if model =~ /[+]PMSF/i
  model.sub!(/[+]PMSF/, '') 
  pmsf = 'relaxed'
end

if ref_treefile.nil?
  STDERR.puts "-r ref_treefile not given! Exiting ......"
  exit(1)
else
  treefile = ref_treefile
end


# Initialize logging
mkdir_with_force(outdir, is_force)
init_logging(outdir)


#########################################
begin
  outdirs[:basics] = File.join(outdir, 'julia')
  outdirs[:bl] = File.join(outdir, 'bl')
  outdirs[:inBV] = File.join(outdir, 'inBV')
  outdirs[phylo_prog] = File.join(outdir, phylo_prog.to_s)
  outdirs.values.each { |sub_outdir| mkdir_with_force(sub_outdir, is_force) }

  branchout_matrix = File.join(outdirs[:bl], 'branch_out.matrix')
  in_BV = File.join(outdir, 'in.BV')
  # parse mfm and regular model: LG+C20 will be parsed as mfm=C20 & non_mfm=LG
  non_mfm, mfm = parse_model_name(model)

  if not iqtree_indir.nil?
    FileUtils.rm_rf(outdirs[:iqtree])
    FileUtils.cp_r(iqtree_indir, outdirs[:iqtree])
  else
    case pmsf
      when 'relaxed'
        iqtree_add_arg = iqtree_add_arg + " -ft #{outdirs[:iqtree]}/iqtree.treefile" # using guide tree, which will generate .sitefreq
        julia_bl_add_arg = "--pmsf #{outdirs[:iqtree]}/iqtree.sitefreq"
      when 'stringent'
        ;
    end
    case phylo_prog
      when :iqtree
        run_iqtree(treefile, model, pmsf, seqfile, outdirs[:iqtree], blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
        out_treefile = File.join(outdirs[:iqtree], 'iqtree.treefile')
      when :phyml
        old_wd = Dir.pwd
        FileUtils.cp_r([treefile, seqfile], outdirs[:phyml])
        Dir.chdir(outdirs[:phyml])
        c = getCorename(seqfile)
        run_phyml(File.basename(treefile), model, pmsf, File.basename(seqfile), outdirs[:phyml], blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
        out_treefile = File.join(outdirs[:phyml], [c,'phy_phyml_tree.txt'].join('.'))
        Dir.chdir(old_wd)
    end
  end

  log("")

  # do_bl_my_try.R to generate julia_outdir (basics)
  cmd = "Rscript #{GEN_BASICS} -s #{seqfile} -t #{out_treefile} --cpu #{cpu} --julia_outdir #{outdirs[:basics]} --force --type #{st}"
  execute_command(cmd, "Generating basic info of the tree and alignment ......")

  # generate_branch_out_mat.sh
  cmd = "bash #{GEN_BRANCH} -t #{out_treefile} --ref_tree #{ref_treefile} --outdir #{outdirs[:bl]} --force"
  execute_command(cmd, "Parsing the branch order ......")
  
  # julia_bl
  cmd = "#{JULIA} -t #{cpu} #{JULIA_BL} --basics_indir #{outdirs[:basics]} -t #{st} --tree #{out_treefile} -b #{branchout_matrix} -m #{non_mfm} --mfm #{mfm} --outdir #{outdirs[:inBV]} --force #{julia_bl_add_arg} --hessian_type #{hessian_type} 2>#{outdir}/error"
  execute_command(cmd, "Calculating Hessian. Takes long time ......")

  phylip = create_phylip(seqfile) if phylip.nil?
  run_mcmctree(outdirs, tree_indir, phylip, clock, bd, bsn) if is_run_mcmctree
  STDOUT.puts "Done!" if $? == 0
rescue => e
  log_error("Fatal error: #{e.message}")
  log_error(e.backtrace.join("\n"))
  exit(1)
ensure
  $log.close if $log
end


