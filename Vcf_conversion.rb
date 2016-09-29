#Methods

def get_line_from_file(path, line)
result = nil
  File.open(path, "r") do |f|
    while line > 0
      line -= 1
      result = f.gets
      result ||= '' 
      result = result.chomp
    end
  end
  return result
end

def write_header(output_path)
  formato = "##fileformat=VCFv4.2\n"
  colonne = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
  File.open(output_path, "w") { |file| file.write(formato + colonne) }
end

def write_line(output_path, chrom, pos, id, ref, alt, qual, filter, info, format, sample)
   File.open(output_path, "a") { |file| file.write(chrom + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + filter + "\t" + info + "\t" + format + "\t" + sample + "\n") }
end

def find_alt(snp_line, ref)

  alternatives = snp_line[9].to_s.split "/"
  alt = "."
  if (alternatives.length == 2)
    alternatives.each do |a|
      if a.to_s != ref
        alt = a.to_s
      end
    end
  else  
    puts "Error: bihallelic field"
    ref = "X"
  end

  return alt
    
end


def create_line(snp_path, current_line, current_father, current_mother, output)
  
  current_snp_line = current_line + 3

  snp_line = get_line_from_file(snp_path, current_snp_line).split "\t"

  chrom = snp_line[1].to_s 

  pos = (snp_line[2].to_i + 1).to_s 

  id = "."

  ref = snp_line[7].to_s

  alt = find_alt(snp_line, ref)
  
  qual = "."

  filter = "PASS"

  info =  "."

  format = "GT" 

  sample = current_father.to_s + "|" + current_mother.to_s
  
  write_line(output, chrom, pos, id, ref, alt, qual, filter, info, format, sample)

end

def write_body(haplotypes, snp, output)

  current_line = 0
  father = get_line_from_file(haplotypes, 1).split('')
  mother = get_line_from_file(haplotypes, 2).split('')
  
  father.each do |i|
    create_line(snp, current_line, father[current_line], mother[current_line], output)
    current_line += 1
  end

end

def create_vcf(haplotypes_path, snp_path, output_path)

  snp_lines = `wc -l "#{snp_path}"`.strip.split(' ')[0].to_i - 3
  
  if (get_line_from_file(haplotypes_path, 1).length == snp_lines)
    write_header(output_path)
    write_body(haplotypes_path, snp_path, output_path)
  else
    puts "Error: non-matching files"
  end

end


#Vcf file creation

create_vcf(haplotypes_path, snp_path, output_path)
