#! /bin/env ruby

# © Sishuo Wang [2024-2025]. All rights reserved.
puts("© Sishuo Wang [2024-2025]. All rights reserved.", "\n")


#########################################
# 2025-04-16


#########################################
require 'tempfile'
require 'fileutils'


#########################################
$VERBOSE = nil
DIR = File.dirname($0)
ADD_SCRIPTS_DIR = File.join(DIR, "additional_scripts")
$VERBOSE = true

JULIA_BL = File.join(DIR, "julia_bl.jl")
GEN_BASICS = File.join(DIR, 'do_bl_my_try.R')
GEN_BRANCH = File.join(ADD_SCRIPTS_DIR, 'generate_branch_out_mat.sh')

IQTREE = 'iqtree'
JULIA = 'julia'
RUBY = 'ruby'

REGULAR_DIR = File.join(DIR, 'substitution_model', 'regular')
MFM_DIR = File.join(DIR, 'substitution_model', 'mfm')

DO_MCMCTREE = File.expand_path("~/lab-tools/dating/do_mcmctree.rb")
BS_INBV = File.expand_path("~/lab-tools/dating/hessian/create_hessian_by_bootstrapping.rb")
MFATOPHY = File.join(ADD_SCRIPTS_DIR, 'MFAtoPHY.pl')


#########################################
require 'getoptlong'
require 'parallel'
require 'colorize'

require_relative 'additional_scripts/Dir'


#########################################
# parse model_name to get mfm if any
def check_model_num(a, model)
  if a.size == 1
    return(a[0])
  elsif a.size >= 2
    raise "Error! More than one mfm is specified in #{model}"
  else
    return('nothing')
  end
end


def parse_model_name(model)
  non_mfms = Array.new
  mfms = Array.new

  model.split('+').each do |m|
    if %w[dat nex].any?{|suffix| File.exist?(File.join(REGULAR_DIR, m + '.' + suffix)) }
      non_mfms << m
    elsif %w[dat nex].any?{|suffix| File.exist?(File.join(MFM_DIR, m + '.' + suffix)) }
      mfms << m
    end
  end

  non_mfm = check_model_num(non_mfms, model)
  non_mfm = non_mfm == 'nothing' ? 'POISSON' : non_mfm
  mfm = check_model_num(mfms, model)

  return([non_mfm, mfm])
end


def run_mcmctree(outdirs, tree_indir, phylip, clock, bd, bsn)
  outdirs[:date] = File.join(File.dirname(outdirs[:inBV]), 'date')
  mkdir_with_force(outdirs[:date])
  inBV_file = File.join(outdirs[:inBV], 'in.BV')
  add_argu = "--inBV #{inBV_file}"
  cmd = "#{RUBY} #{DO_MCMCTREE} --outdir #{outdirs[:date]} --tree_indir #{tree_indir} -i #{phylip} --prot --clock #{clock} --bd #{bd} --bsn #{bsn} --print 2 #{add_argu} --force --one_codeml"
  p cmd
  ` #{cmd} `
end


def create_phylip(seqfile)
  tmp_dir = Dir.mktmpdir
  `#{MFATOPHY} #{seqfile}; mv #{seqfile}.phy #{tmp_dir}`
  phylip = File.join(tmp_dir, File.basename(seqfile)+'.phy')
  return(phylip)
end


def run_iqtree(treefile, model, pmsf, seqfile, outdirs, blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
  if ! pmsf.nil?
    `#{IQTREE} -te #{treefile} -m #{model} -s #{seqfile} -pre #{outdirs[:iqtree]}/iqtree -redo -wsl -quiet #{iqtree_add_arg0} -T #{cpu}`
  end

  `#{IQTREE} -blfix #{treefile} -m POISSON -s #{seqfile} -pre #{outdirs[:iqtree]}/iqtree -redo -wsl -quiet #{iqtree_add_arg0} -T #{cpu}`
  `#{IQTREE} -te #{outdirs[:iqtree]}/iqtree.treefile -m #{model} -s #{seqfile} -pre #{outdirs[:iqtree]}/iqtree -redo -wsl -quiet #{iqtree_add_arg} -blmin #{blmin} -T #{cpu}`
end


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
iqtree_add_arg0 = ['-mwopt', '-keep-ident', "-mdef #{MFM_DIR}/C2.nex"].join(' ')
iqtree_add_arg = iqtree_add_arg0
julia_bl_add_arg = nil

hessian_type = 'SKT2004'
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
      phylo_prog = value
    when '--phyml'
      phylo_prog = :phyml
    when '--iqtree'
      phylo_prog = :iqtree
    when '--iqtree_indir'
      iqtree_indir = value
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


#########################################
mkdir_with_force(outdir, is_force)

outdirs[:basics] = File.join(outdir, 'julia')
outdirs[:iqtree] = File.join(outdir, 'iqtree')
outdirs[:bl] = File.join(outdir, 'bl')
outdirs[:inBV] = File.join(outdir, 'inBV')
outdirs.values.map{|sub_outdir| mkdir_with_force(sub_outdir, is_force) }


branchout_matrix = File.join(outdirs[:bl], 'branch_out.matrix')
in_BV = File.join(outdir, 'in.BV')

# parse mfm and regular model: LG+C20 will be parsed as mfm=C20 & non_mfm=LG
non_mfm, mfm = parse_model_name(model)


if not iqtree_indir.nil?
  FileUtils.rm_rf(outdirs[:iqtree])
  FileUtils.cp_r(iqtree_indir, outdirs[:iqtree])
  #`cp -r #{iqtree_indir}/* #{outdirs[:iqtree]}`
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
      run_iqtree(treefile, model, pmsf, seqfile, outdirs, blmin, iqtree_add_arg0, iqtree_add_arg, cpu)
    when :phyml
      run_phyml()
  end

end

out_treefile = File.join(outdirs[:iqtree], 'iqtree.treefile')


#########################################
# do_bl_my_try.R to generate julia_outdir (basics)
puts "Rscript #{GEN_BASICS} -s #{seqfile} -t #{out_treefile} --cpu #{cpu} --julia_outdir #{outdirs[:basics]} --force --type #{st}"
`Rscript #{GEN_BASICS} -s #{seqfile} -t #{out_treefile} --cpu #{cpu} --julia_outdir #{outdirs[:basics]} --force --type #{st}`

# generate_branch_out_mat.sh
puts "bash #{GEN_BRANCH} -t #{out_treefile} --ref_tree #{ref_treefile} --outdir #{outdirs[:bl]} --force"
`bash #{GEN_BRANCH} -t #{out_treefile} --ref_tree #{ref_treefile} --outdir #{outdirs[:bl]} --force`

# julia_bl
cmd = "#{JULIA} -t #{cpu} #{JULIA_BL} --basics_indir #{outdirs[:basics]} -t #{st} --tree #{out_treefile} -b #{branchout_matrix} -m #{non_mfm} --mfm #{mfm} --outdir #{outdirs[:inBV]} --force #{julia_bl_add_arg} --hessian_type #{hessian_type} 2>#{outdir}/error"
puts cmd
`#{cmd}`


if phylip.nil?
  phylip = create_phylip(seqfile)
end

run_mcmctree(outdirs, tree_indir, phylip, clock, bd, bsn) if is_run_mcmctree


