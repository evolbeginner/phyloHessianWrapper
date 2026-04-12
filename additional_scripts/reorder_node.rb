#! /bin/env ruby


###############################################################
require 'getoptlong'
require 'bio-nwk'
require 'Dir'


###############################################################
def get_ordered_tips(infile)
  ordered_tips = Array.new
  tree = nil

  in_fh = File.open(infile, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line.scan(/([^(),: ]+):/) do |i|
      ordered_tips << i[0]
    end
    tree = getTreeObjFromNwkString(line)
    break
  end
  in_fh.close
  #ordered_tips.sort_by!(&:downcase)

  return([ordered_tips, tree])
end


def get_child(tree, node, bls, name2bl, name2order, order2names, name2edge, is_name2order=false)
  if node != tree.root
    bls << tree.distance(node, tree.parent(node))
    #edge = tree.get_edge(node, tree.parent(node))
    #edge.distance = 999
    tips = tree.tips(node)

    tip_names = tips.map{|i|i.name.gsub(' ', '_')}.sort_by(&:downcase)
    complement_tip_names = (ORDERED_TIP_NAMES - tree.tips(node).map{|t|t.name.gsub(' ', '_')}).sort_by(&:downcase)
    #p tree.tips(node).map{|t|t.name}

    name2bl[tip_names] = bls[-1]
    name2bl[complement_tip_names]= bls[-1]

    edge = tree.get_edge(node, tree.parent(node))
    name2edge[tip_names] = edge
    name2edge[complement_tip_names]= edge

    if is_name2order
      order = order2names.size # the 1st is zero
      name2order[tip_names] = order
      name2order[complement_tip_names] = order
      order2names[order] = [tip_names, complement_tip_names]
    end
  end


  if node.isTip?(tree)
    ;
  else
    children = tree.children(node)
    ordered_children = children.sort_by{|child| tree.tips(child).map{|i|ORDERED_TIP_NAMES.index(i.name.gsub(' ', '_'))}.min}
    ordered_children.each do |i|
      get_child(tree, i, bls, name2bl, name2order, order2names, name2edge, is_name2order)
    end
  end
end


###############################################################
def eh(x)
  eh0 = 1e-8
  #(eh0 * (x.abs + 1)) ** 0.67
  eh0 * (x.abs)
end


###############################################################
infile = nil
ref_tree_file = nil
is_output_branch = false
bl_outdir = nil
is_force = false

name2bl = Hash.new
name2edge = Hash.new
name2order, order2names = [Hash.new, Hash.new]


###############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
  ['--ref', GetoptLong::REQUIRED_ARGUMENT],
  ['--output_branch', GetoptLong::NO_ARGUMENT],
  ['--bl_outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
    when '--ref'
      ref_tree_file = value
    when '--output_branch'
      is_output_branch = true
    when '--bl_outdir'
      bl_outdir = value
    when '--force'
      is_force = true
  end
end

if ! bl_outdir.nil?
  mkdir_with_force(bl_outdir, is_force)
end


###############################################################
ORDERED_TIP_NAMES, ref_tree = get_ordered_tips(ref_tree_file)

bls = Array.new
ref_tree = get_child(ref_tree, ref_tree.root, bls, name2bl, name2order, order2names, name2edge, true)
NAME2ORDER = name2order
ORDER2NAMES = order2names

name2bl = Hash.new


###############################################################
trees = getTreeObjs(infile)

trees.each do |tree|
  bls = Array.new
  root = tree.root
  get_child(tree, root, bls, name2bl, name2order, order2names, name2edge, false)

  bls2 = Array.new
  ORDER2NAMES.each_pair do |order, names|
    names = order2names[order]
    bls2 << names.map{|name| name2bl[name] }.compact[0]

    if not bl_outdir.nil?
      count = 0
      name = names[0]
      old_dist = name2edge[name].distance
      #change_bl(name, name2edge)
      sub_outdir = File.join(bl_outdir, order.to_s)
      mkdir_with_force(sub_outdir, is_force)

      [[1, :+], [1, :-]].each do |delta, op|
        count += 1
        sub_sub_outdir = File.join(sub_outdir, count.to_s)
        mkdir_with_force(sub_sub_outdir)
        #factor = old_dist > 1e-6 ? 1e-6 : old_dist*0.1
        step = eh(old_dist)
        name2edge[name].distance = old_dist.send(op, step)
        File.write(File.join(sub_sub_outdir, 'bl_changed.treefile'), tree.cleanNewick + "\n")
        File.write(File.join(sub_sub_outdir, 'step_size'), step.to_s + "\n")
      end

      name2edge[name].distance = old_dist
    end

    if is_output_branch
      puts [[names.map{|a|a.join('-')}].join(','), bls2[-1]].join("\t")
    else
      if names.map{|name| name2bl[name] }.any?{|i|i.nil?}
        p names
        p names.map{|name| name2bl[name] }
      end
    end
    #p [bls[-1], names]
  end
  puts bls2.join(' ') if not is_output_branch
end


