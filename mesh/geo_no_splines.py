from sys import argv, exit
import re

if len(argv<2):
    print(f"Usage: {argv[0]} <entrada.geo>")
    exit(1)

geo_filename = sys.argv[1]

line_count = 1
loop_count = 1
loop_list = []
skip_next_line = False

# List of loops to be deleted (islands, etc)
skip_loops = [2, 3, 4, 5, 7, 8, 21]

def get_spline_id(line):
    string="Spline(IL+aa) = "
    expression="\w*\((IL\+aa)\)"
    result = re.match(expression, string)
    assert(result)
    # print("OK")
    id = result.group(1)
    return id

def process_line(line):
    global line_count
    global loop_count
    global loop_list

    if not loop_count in skip_loops:

        pts = [ s[0:-1] for s in line.split() ] # Dividir y quitar coma final
        loop_id = get_spline_id(line)
        pts = pts[2:] # Eliminar título y " = "
        pts # Ahí estarán todos los identificadores de los puntos
        pts[0] = pts[0][1:] # Quitar "{" al primero
        pts[-1] = pts[-1][:-1] # y "}" al último
        # print(f"Puntos: {pts}")
        n=len(pts)
        first_line = line_count
        line_ids = []
        for i in range(1,n):
            p0 = pts[i-1]
            p1 = pts[i]
            line_ids.append(line_count)
            print(f"Line(" + str(line_count) + ") = {" + str(p0) + ", " + str(p1) + "};")
            line_count = line_count + 1
        last_line = line_count-1
        lines = f"{line_ids[0]}"
        for l_id in line_ids[1:]:
            lines = f"{lines}, {l_id}"
        print(f"Line Loop({loop_count}) = {{ {lines} }};")
        loop_list.append(loop_count)

    loop_count = loop_count + 1

if __name__ == "__main__":
    for line in open(geo_filename, "r"):
        if line[:6] == "Spline":
            process_line(line)
            skip_next_line = True
        elif line[:13] == "Plane Surface":
            loops = f"{loop_list[0]}"
            for l_id in loop_list[1:]:
                loops = f"{loops}, {l_id}"
            print( f"Plane Surface(IS) = {{{loops}}};");
        else:
            if skip_next_line:
                skip_next_line = False
            else:
                print(line, end="")



# line = "Spline(IL+0) = {IP+0, IP+1, IP+2, IP+3, IP+4, IP+5, IP+6, IP+7, IP+8, IP+9, IP+10, IP+11, IP+12, IP+13, IP+14, IP+15, IP+16, IP+17, IP+18, IP+19, IP+20, IP+21, IP+22, IP+23, IP+24, IP+25, IP+26, IP+27, IP+28, IP+29, IP+30, IP+31, IP+32, IP+33, IP+34, IP+35, IP+36, IP+37, IP+38, IP+39, IP+40, IP+41, IP+42, IP+43, IP+44, IP+45, IP+46, IP+47, IP+48, IP+49, IP+50, IP+51, IP+52, IP+53, IP+54, IP+55, IP+56, IP+57, IP+58, IP+59, IP+60, IP+61, IP+62, IP+63, IP+64, IP+65, IP+66, IP+67, IP+68, IP+69, IP+70, IP+71, IP+72, IP+73, IP+74, IP+75, IP+76, IP+77, IP+78, IP+79, IP+80, IP+81, IP+82, IP+83, IP+84, IP+85, IP+86, IP+87, IP+88, IP+89, IP+90, IP+91, IP+92, IP+93, IP+94, IP+95, IP+96, IP+97, IP+98, IP+99, IP+100, IP+101, IP+102, IP+103, IP+104, IP+105, IP+106, IP+107, IP+108, IP+109, IP+110, IP+111, IP+112, IP+113, IP+114, IP+115, IP+116, IP+117, IP+118, IP+119, IP+120, IP+121, IP+122, IP+123, IP+124, IP+125, IP+126, IP+127, IP+128, IP+129, IP+130, IP+131, IP+132, IP+133, IP+134, IP+135, IP+136, IP+137, IP+138, IP+139, IP+140, IP+141, IP+142, IP+143, IP+144, IP+145, IP+146, IP+147, IP+148, IP+149, IP+150, IP+151, IP+152, IP+153, IP+154, IP+155, IP+156, IP+157, IP+158, IP+159, IP+160, IP+161, IP+162, IP+163, IP+164, IP+165, IP+166, IP+167, IP+168, IP+169, IP+170, IP+171, IP+172, IP+173, IP+174, IP+175, IP+176, IP+177, IP+0};"
