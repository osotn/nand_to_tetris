#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include <SDL.h>
#include <SDL2/SDL_ttf.h>

#pragma GCC diagnostic ignored "-Wwrite-strings"

typedef uint8_t bit_t;

typedef struct {
	bit_t bit[16];
} bit16_t;

const bit16_t bit16_zeros = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
const bit16_t bit16_ones  = {{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
const bit16_t bit16_one   = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};

const char* bit_to_string(bit_t bit)
{
	return bit ? "1" : "0";
}

void print_bit(char * pstr, bit_t bit)
{
        printf("%s = %s\n", pstr, bit_to_string(bit));
}

void int16_to_bit16(int16_t val, bit16_t* bit16)
{
	uint16_t u = (int16_t)val;
	int i;

	for (i=0; i<16; i++) {
		bit16->bit[i] = (u & 0x1) ? 1 : 0;
	      	u >>= 1;
	}	
}

int16_t bit16_to_int16(bit16_t* bit16)
{
	uint16_t u = 0;
	int i;

	for (i=0; i<16; i++) {
		u <<= 1;
		u |= bit16->bit[15-i] ? 0x1 : 0x0;
	}

	return (int16_t)u;
}

#define BIT16_SIZE 19

char* bit16_to_string(char* str, bit16_t* pbit16)
{
	int i;
	str[0] = 'b';
	str[1] = '\'';
	for (i=0; i<16; i++)
		str[2+i] = pbit16->bit[15-i] ? '1' : '0';
	str[18] = '\0';
	return str;
}

int16_t string_bit16_to_int16(char* str)
{
    int i;
    uint16_t u = 0;
    for (i=0; i<16; i++) {
        u <<= 1;
        u |= (str[2+i] != '0') ? 0x1 : 0x0;
    }

    return (int16_t)u;
}

char* int16_to_string(char* str, int16_t v)
{
    bit16_t bit16;
    int16_to_bit16(v, &bit16);
    return bit16_to_string(str, &bit16);
}

char* uint16_bits_to_string(char* str, uint16_t u16, int n)
{
    int i;
    str[0] = 'b';
    str[1] = '\'';
    bit16_t bit16;
    int16_to_bit16((int16_t)u16, &bit16);

    for (i = 0; i < n; i++)
        str[2+i] = bit16.bit[n-1-i] ? '1' : '0';
    str[2+n] = '\0';
    return str;
}

void print_bit16(char * pstr, bit16_t* pbit16)
{
	char str[BIT16_SIZE];
	printf("%s = %s\n", pstr, bit16_to_string(str, pbit16));
}

typedef struct {
	bit_t  bit;
} dff_t;

uint32_t num_dff = 0;
void dff_init(dff_t* pdff)
{
	pdff->bit = 0;
	num_dff++;
}


bit_t dff_get_output(dff_t* pdff)
{
	return pdff->bit;
}

void print_dff(dff_t* pdff)
{
	print_bit("dff.bit", pdff->bit);
}

void dff_run(dff_t *pdff, bit_t in_bit)
{
	pdff->bit = in_bit;
}

typedef struct {
	dff_t dff[16];
} dff16_t;

void dff16_init(dff16_t* pdff16) {
	int i;
	for (i = 0; i < 16; i++)
		dff_init(&(pdff16->dff[i]));
}

bit16_t dff16_get_output(dff16_t* pdff16) {
	bit16_t bit16;
	int i;

	for (i=0; i<16; i++)
        	bit16.bit[i] = dff_get_output(&(pdff16->dff[i]));

	return bit16;
}

void print_dff16(dff16_t* pdff16)
{
	bit16_t bit16 = dff16_get_output(pdff16);
	print_bit16("dff16.bit16", &bit16);
}

void dff16_run(dff16_t* pdff16, bit16_t* pbit16)
{
	int i;

	for (i=0; i<16; i++)
		dff_run(&(pdff16->dff[i]), pbit16->bit[i]);	
}	

uint32_t num_nand = 0;
void nand(bit_t a, bit_t b, bit_t* out)
{
	a = !!a;
	b = !!b;
	*out = !(a && b);
	num_nand++;
}

uint32_t num_and = 0;
void gate_and(bit_t a, bit_t b, bit_t* out)
{
	a = !!a;
	b = !!b;
	*out = a && b;
	num_and++;
}

uint32_t num_or = 0;
void gate_or(bit_t a, bit_t b, bit_t* out)
{
	a = !!a;
	b = !!b;
	*out = a || b;
	num_or++;
}

uint32_t num_not = 0;
bit_t gate_not(bit_t in, bit_t* out)
{
	*out = ! in;
	num_not++;
	return *out;
}

void not16(bit16_t* in16, bit16_t* out16)
{
	int i;
	for (i=0; i<16; i++)
		gate_not(in16->bit[i], &(out16->bit[i]));
}


void and16(bit16_t* a16, bit16_t* b16, bit16_t* out16)
{
	int i;
	for (i=0; i<16; i++)
		gate_and(a16->bit[i], b16->bit[i], &(out16->bit[i]));
}

void or16_out1(bit16_t* in16, bit_t* out)
{
	*out = 0;
	int i;
	for (i=0; i<16; i++)
		gate_or(in16->bit[i], *out, out);
}

void or3_out1(bit_t a0, bit_t a1, bit_t a2, bit_t* out)
{
	*out = 0;
	gate_or(*out, a0, out);
	gate_or(*out, a1, out);
	gate_or(*out, a2, out);
}

void gate_xor(bit_t a, bit_t b, bit_t* out)
{
	bit_t na, nb, a_and_nb, b_and_na;

	gate_not(a, &na);
	gate_not(b, &nb);
	gate_and(a, nb, &a_and_nb);
	gate_and(b, na, &b_and_na);
	gate_or(a_and_nb, b_and_na, out);
}

void mux(bit_t a, bit_t b, bit_t sel, bit_t* out)
{
	bit_t nsel, a_and_nsel, b_and_sel;

	gate_not(sel, &nsel);
	gate_and(a, nsel, &a_and_nsel);
	gate_and(b, sel, &b_and_sel);
	gate_or(a_and_nsel, b_and_sel, out);
}

void dmux(bit_t in, bit_t sel, bit_t* a, bit_t*b)
{
    bit_t nsel;

	gate_not(sel, &nsel);
	gate_and(in, nsel, a);
	gate_and(in, sel, b);
}

void dmux4way(bit_t in, bit_t sel0, bit_t sel1, bit_t* a, bit_t*b, bit_t*c, bit_t*d)
{
    bit_t ab, cd;

	dmux(in, sel1, &ab, &cd);
	dmux(ab, sel0, a, b);
	dmux(cd, sel0, c, d);

}


void mux16(bit16_t* a, bit16_t* b, bit_t sel, bit16_t* out)
{
	int i;
	for (i=0; i<16; i++)
		mux(a->bit[i], b->bit[i], sel, &(out->bit[i]));
}

void mux4way16(bit16_t* a, bit16_t* b, bit16_t* c, bit16_t* d, bit_t sel0, bit_t sel1, bit16_t* out)
{
    bit16_t ab, cd;

    mux16(a, b, sel0, &ab);
    mux16(c, d, sel0, &cd);
    mux16(&ab, &cd, sel1, out);
}

void hadd(bit_t a, bit_t b, bit_t* sum, bit_t* carry)
{
	gate_xor(a, b, sum);
	gate_and(a, b, carry);
}

void add(bit_t a, bit_t b, bit_t in_carry, bit_t* sum, bit_t* out_carry)
{
	bit_t sum1, carry1, carry2;

	hadd(a, b, &sum1, &carry1);
	hadd(in_carry, sum1, sum, &carry2);
	gate_or(carry1, carry2, out_carry);
}

void add16(bit16_t* a, bit16_t* b, bit16_t* sum)
{
	bit_t carry = 0;

	int i;
	for (i=0; i<16; i++)
		add(a->bit[i], b->bit[i], carry, &(sum->bit[i]), &carry);
}

void alu(bit16_t* x16, bit16_t* y16, bit_t zx, bit_t nx, bit_t zy, bit_t ny, bit_t f, bit_t nf, bit16_t* out, bit_t* nz, bit_t* ng)
{
	bit16_t zero = bit16_zeros;
	bit16_t x16_1, x16_2, x16_3;
	bit16_t y16_1, y16_2, y16_3;
	bit16_t out16_1, out16_2, out16_3, out16_4;

/*
        zx | nx | zy | ny | f | no | out
         1 |  0 |  1 |  0 | 1 |  0 |  0           out = 0 + 0
         1 |  1 |  1 |  1 | 1 |  1 |  1           out = !(-1 + -1) = !(-2) = 1
         1 |  1 |  1 |  0 | 1 |  0 |  -1          out = -1 + 0 = -1
         0 |  0 |  1 |  1 | 0 |  0 |  x           out = x & -1 = x
         1 |  1 |  0 |  0 | 0 |  0 |  y           out =  -1 & y = y
         0 |  0 |  1 |  1 | 0 |  1 |  !x
         1 |  1 |  0 |  0 | 0 |  1 |  !y
         0 |  0 |  1 |  1 | 1 |  1 |  -x
         1 |  1 |  0 |  0 | 1 |  1 |  -y
         0 |  1 |  1 |  1 | 1 |  1 |  x+1
         1 |  1 |  0 |  1 | 1 |  1 |  y+1
*/


	mux16(x16, &zero, zx, &x16_1);
	not16(&x16_1, &x16_2);
	mux16(&x16_1, &x16_2, nx, &x16_3);
	
	mux16(y16, &zero, zy, &y16_1);
	not16(&y16_1, &y16_2);
	mux16(&y16_1, &y16_2, ny, &y16_3);

	and16(&x16_3, &y16_3, &out16_1);
	add16(&x16_3, &y16_3, &out16_2);

	mux16(&out16_1, &out16_2, f, &out16_3);

	not16(&out16_3, &out16_4);
	mux16(&out16_3, &out16_4, nf, out);

    or16_out1(out, nz);
	*ng = out->bit[15];
}

typedef struct {
	dff16_t dff16;
	char* name;
} reg16_t;

void reg16_init(reg16_t* reg16, char* name) {
	dff16_init(&(reg16->dff16));
	reg16->name = name;
}

bit16_t reg16_get_output(reg16_t* reg16) {
	bit16_t bit16 = dff16_get_output(&(reg16->dff16));
	return bit16;
}

void print_reg16(reg16_t* reg16)
{
	bit16_t bit16 = reg16_get_output(reg16);
	char str[BIT16_SIZE];
	printf("%s:reg16.bit16 = %s\n", reg16->name, bit16_to_string(str, &bit16));
}

void reg16_run(reg16_t* reg16, bit16_t* in16, bit_t load)
{
	bit16_t bit16;
        bit16_t bit16_z1 = dff16_get_output(&(reg16->dff16));
	mux16(&bit16_z1, in16, load, &bit16);
	dff16_run(&(reg16->dff16), &bit16);
}	

typedef struct {
	dff16_t dff16;
} pc16_t;

void pc16_init(pc16_t* pc16) {
	dff16_init(&(pc16->dff16));
}

bit16_t pc16_get_output(pc16_t* pc16) {
	bit16_t bit16 = dff16_get_output(&(pc16->dff16));
	return bit16;
}

void print_pc16(pc16_t* pc16)
{
	bit16_t bit16 = pc16_get_output(pc16);
	print_bit16("pc16.bit16", &bit16);
}

void pc16_run(pc16_t* pc16, bit_t reset, bit16_t* in16, bit_t load)
{
	bit16_t bit16_1, bit16_2;
        bit16_t bit16_z1 = dff16_get_output(&(pc16->dff16));
	bit16_t one = bit16_one;
	bit16_t zero = bit16_zeros;
	bit16_t inc16;

	add16(&bit16_z1, &one, &inc16);
	mux16(&inc16, in16, load, &bit16_1);
	mux16(&bit16_1, &zero, reset, &bit16_2);
	dff16_run(&(pc16->dff16), &bit16_2);
}	

#define MEM_16K 0x4000

typedef struct {
    int16_t mem[MEM_16K];
} ram16k_t;

void ram16k_init(ram16k_t* ram16k)
{
    int i;
    for (i=0; i<MEM_16K; i++)
        ram16k->mem[i] = 0;
}

bit16_t ram16k_get_output(ram16k_t* ram16k, bit16_t* addr)
{
    bit16_t out = bit16_zeros;
    uint16_t mem_addr = (uint16_t)bit16_to_int16(addr);
    if (mem_addr < MEM_16K)
        int16_to_bit16(ram16k->mem[mem_addr], &out);
    return out;
}

void ram16k_run(ram16k_t* ram16k, bit16_t* in, bit_t load, bit16_t* addr)
{
    uint16_t mem_addr = (uint16_t)bit16_to_int16(addr);

    if (mem_addr < MEM_16K) {
        if (load)
            ram16k->mem[mem_addr] = bit16_to_int16(in);
    }
}

void ram16k_dump(ram16k_t* ram16k, uint16_t addr, int n)
{
    printf("RAM16k dump, addr = 0x%x, n = %d:\n", addr, n);
    int i;
    char str[BIT16_SIZE];
    for (i=0; i<n; i++, addr++) {
        if (addr < MEM_16K) {
            int16_t v = ram16k->mem[addr];
            printf(" [0x%d] = %d, 0x%x, %s\n", addr, v, v, int16_to_string(str, v));
        }
    }
}

#define MEM_8K 0x2000

typedef struct {
    int16_t mem[MEM_8K];
} screen_t;

void screen_init(screen_t* screen)
{
    int i;
    for (i=0; i<MEM_8K; i++)
        screen->mem[i] = 0;
}

bit16_t screen_get_output(screen_t* screen, bit16_t* addr)
{
    bit16_t out = bit16_zeros;
    uint16_t mem_addr = (uint16_t)bit16_to_int16(addr);
    if (mem_addr < MEM_8K)
        int16_to_bit16(screen->mem[mem_addr], &out);
    return out;
}

uint16_t screen_debug = 0;

void screen_run(screen_t* screen, bit16_t* in, bit_t load, bit16_t* addr)
{
    uint16_t mem_addr = (uint16_t)bit16_to_int16(addr);

    if (screen_debug) {
                char str[BIT16_SIZE];
                printf("screen debug: mem_addr = %d, <== %s\n",
                        mem_addr, bit16_to_string(str, addr));
    }

    if (mem_addr < MEM_8K) {
        if (load) {
            screen->mem[mem_addr] = bit16_to_int16(in);
            if (screen_debug) {
                char str[BIT16_SIZE];
                printf("screen debug: Load mem_addr = %d, in = %s => [%d]=0x%x\n",
                        mem_addr, bit16_to_string(str, in), mem_addr, screen->mem[mem_addr]);
            }
        }
    }
}

#define SCREEN_ROW_NUM 256
#define SCREEN_COLUMN_NUM 512
#define SCREEN_PIXEL_BITS 1
#define SCREEN_COLUMN_WORDS (512 * SCREEN_PIXEL_BITS / 16)

uint16_t screen_get_pixel(screen_t* screen, int r, int c)
{
    if (r < 0 || r >= SCREEN_ROW_NUM ||
        c < 0 || c >= SCREEN_COLUMN_NUM) {
        return 0;
    }

    int addr = r * SCREEN_COLUMN_WORDS + c / 16;
    if (addr >= MEM_8K) {
        return 0;
    }
    uint16_t word = screen->mem[addr];

    return  (word & (1u << ((uint16_t)c % 16))) ? 1 : 0;
}

void screen_dump(screen_t* screen, int r_i, int r_n, int c_i, int c_n)
{
    printf("SCREEN dump, Row i=%d n=%d; Column i=%d n=%d\n", r_i, r_n, c_i, c_n);
    int r, c;
    for (r=0; r < r_n; r++) {
        for (c=0; c < c_n; c++) {
            printf("%s", screen_get_pixel(screen, r, c) ? "*" : " ");
        }
        printf("\n");
    }
}

int screen_save_bmp(screen_t* screen, char* file_name)
{
    struct {
        uint8_t     signature[2];
        uint32_t    file_size;
        uint32_t    reserved_0;
        uint32_t    data_offset;
    } __attribute__((packed)) bmp_header;
    struct {
        uint32_t     size;
        uint32_t     width;
        uint32_t     height;
        uint16_t     planes;
        uint16_t     bit_per_pixel;
        uint32_t     compression;
        uint32_t     image_size;
        uint32_t     x_pixels_per_m;
        uint32_t     y_pixels_per_m;
        uint32_t     colors_used;
        uint32_t     colors_important;
    } __attribute__((packed)) bmp_info_header;

    uint8_t rgb_black[3] = {0x00, 0x00, 0x00};
    uint8_t rgb_white[3] = {0xFF, 0xFF, 0xFF};

    FILE* file = fopen(file_name, "w");

    if (file == NULL) {
        printf("screen save bmp: - err to create %s file\n", file_name);
        return -1;
    }

    memset(&bmp_header, 0, sizeof(bmp_header));
    bmp_header.signature[0] = 'B';
    bmp_header.signature[1] = 'M';
    bmp_header.data_offset = sizeof(bmp_header) + sizeof(bmp_info_header);

    memset(&bmp_info_header, 0, sizeof(bmp_info_header));
    bmp_info_header.size = sizeof(bmp_info_header);
    bmp_info_header.width = SCREEN_COLUMN_NUM;
    bmp_info_header.height = SCREEN_ROW_NUM;
    bmp_info_header.planes = 1;
    bmp_info_header.bit_per_pixel = 24;
    bmp_info_header.compression = 0;

    fwrite(&bmp_header, sizeof(bmp_header), 1, file);
    fwrite(&bmp_info_header, sizeof(bmp_info_header), 1, file);

    int r, c;
    for (r= 0; r < SCREEN_ROW_NUM; r++) {
        for (c=0; c < SCREEN_COLUMN_NUM; c++) {
            uint8_t* pColor = screen_get_pixel(screen, SCREEN_ROW_NUM-1-r, c) ? 
                                    rgb_black : rgb_white;
            fwrite(pColor, 1, 3, file);
        }
    }

    fclose(file);
    return 0;
}

typedef struct {
    uint16_t code;
} keyboard_t;

void keyboard_init(keyboard_t* keyboard)
{
    keyboard->code = 0;
}

bit16_t keyboard_get_output(keyboard_t* keyboard)
{
    bit16_t out;
    int16_to_bit16((int16_t)keyboard->code, &out);
    return out;
}

void keyboard_set_code(keyboard_t* keyboard, uint16_t code)
{
    keyboard->code = code;
}

void keyboard_run(keyboard_t* keyboard)
{
    // Nothing
}

uint16_t memory_debug = 0;

typedef struct {
    ram16k_t   ram16k;
    screen_t   screen;
    keyboard_t keyboard;
} memory_t;

void memory_init(memory_t* memory)
{
    ram16k_init(&(memory->ram16k));
    screen_init(&(memory->screen));
    keyboard_init(&(memory->keyboard));
}

bit16_t memory_get_output(memory_t* memory, bit16_t* addr)
{
    bit16_t out;
    bit_t sel0 = addr->bit[13];
    bit_t sel1 = addr->bit[14];
    bit16_t out_ram16k = ram16k_get_output(&(memory->ram16k), addr);
    bit16_t out_screen = screen_get_output(&(memory->screen), addr);
    bit16_t out_keyboard = keyboard_get_output(&(memory->keyboard));
    mux4way16(&out_ram16k, &out_ram16k, &out_screen, &out_keyboard, sel0, sel1, &out);
    return out;
}

void memory_run(memory_t* memory, bit16_t* in, bit_t load, bit16_t* addr)
{
    bit_t sel0 = addr->bit[13];
    bit_t sel1 = addr->bit[14];
    bit_t a, b, c, d;
    bit_t ab;

    dmux4way(load, sel0, sel1, &a, &b, &c, &d);
    gate_or(a, b, &ab);

    if (memory_debug) {
        char str1[BIT16_SIZE];
        char str2[BIT16_SIZE];
        printf("mem debug: addr = %s, in = %s, load = %s; ab = %s, c = %s, d = %s\n",
               bit16_to_string(str1, addr), bit16_to_string(str2, in), bit_to_string(load),
               bit_to_string(ab), bit_to_string(c), bit_to_string(d));
    }

    bit16_t addr14 = *addr;
    addr14.bit[15] = 0;
    addr14.bit[14] = 0;
    ram16k_run(&(memory->ram16k), in, ab, &addr14);
    addr14.bit[13] = 0;
    screen_run(&(memory->screen), in, c, &addr14);
    keyboard_run(&(memory->keyboard));
}

#define MEM_32K 0x8000

typedef struct {
    int16_t mem[MEM_32K];
} rom32k_t;

void rom32k_init(rom32k_t* rom32k)
{
    int i;
    for (i = 0; i < MEM_32K; i++)
        rom32k->mem[i] = 0;
}

bit16_t rom32k_get_output(rom32k_t* rom32k, bit16_t* addr)
{
    uint16_t mem_addr = (uint16_t)bit16_to_int16(addr);
    bit16_t out = bit16_zeros;

    if (mem_addr < MEM_32K)
        int16_to_bit16(rom32k->mem[mem_addr], &out);

    return out;
}

void rom32k_load_from_buf(rom32k_t* rom32k, char**pbuf, int n)
{
    int i;
    for (i = 0; i < n; i++)
        rom32k->mem[i] = string_bit16_to_int16(pbuf[i]);
}

int rom32k_load(rom32k_t *rom32k, char* file_name)
{
    FILE* file = fopen(file_name, "r");

    if (file == NULL) {
        printf("rom32k: load - error to open file %s\n", file_name);
        return -1;
    }

    char line[16+1+1];
    int  l = 0;
    int  len;
    while (fgets(line, sizeof(line), file) != NULL) {
        if ((len = strlen(line)) != (sizeof(line)-1)) {
            printf("rom32k: load file %s - error line %d, len = %d != 16+1\n",
                    file_name, l, len);
            fclose(file);
            return -1;
        }
        uint16_t u = 0;
        int i;
        for (i=0; i<16; i++) {
            u <<= 1;
            u |= (line[i] != '0') ? 0x1 : 0x0;
        }
        if (l >= MEM_32K) {
            printf("rom32k: load file %s - warning line number more than %d\n",
                    file_name, MEM_32K);
            fclose(file);
            return -1;
        }
        rom32k->mem[l++] = (int16_t)u;
    }
    printf("rom32k: load %s, %d instructions\n", file_name, l);
    fclose(file);
    return 0;
}

typedef struct {
	reg16_t reg_a16, reg_d16;
	pc16_t  pc16;
	// reference feedback data
	memory_t* ext_memory;
	rom32k_t* ext_rom;
} cpu_t;

void cpu_init(cpu_t* cpu, memory_t* ext_memory, rom32k_t* ext_rom)
{
    cpu->ext_memory = ext_memory;
    cpu->ext_rom    = ext_rom;

	reg16_init(&(cpu->reg_a16), "cpu:reg_a16");
	reg16_init(&(cpu->reg_d16), "cpu:reg_d16");
	pc16_init(&(cpu->pc16));
}

#define INSTRUCTION_SIZE 80

#define AR_SIZEOF(ar) (sizeof(ar)/sizeof(ar[0])) 

/* C instruction */
char* asm_a_dest_name[] = {
        "NULL",     // 0 0 0
        "M",        // 0 0 1
        "D",        // 0 1 0
        "DM",       // 0 1 1
        "A",        // 1 0 0
        "AM",       // 1 0 1
        "AD",       // 1 1 0
        "ADM"       // 1 1 1
};

struct {char* name; uint8_t bits;} asm_a_dest_st[] = {
        { "NULL", 0x0  },
        { "A",    0x4  },
        { "D",    0x2  },
        { "M",    0x1  },
        { "AD",   0x6  }, { "DA",   0x6  },
        { "AM",   0x5  }, { "MA",   0x5  },
        { "DM",   0x3  }, { "MD",   0x3  },
        { "ADM",  0x7  }, { "AMD",  0x7  }, { "DAM",  0x7  }, { "DMA",  0x7  }, { "MAD",  0x7  }, { "MDA",  0x7  }
};

char* asm_a_jump_name[] = {
        "NULL",     // 0 0 0
        "JGT",      // 0 0 1
        "JEQ",      // 0 1 0
        "JGE",      // 0 1 1
        "JLT",      // 1 0 0
        "JNE",      // 1 0 1
        "JLE",      // 1 1 0
        "JMP"       // 1 1 1
};
struct {char* name; uint8_t bits;} asm_a_comp_st[] = {
        { "D&A", 0x0  }, { "D&M", 0x40 }, 
        { "D+A", 0x02 }, { "D+M", 0x42 },
        { "A-D", 0x07 }, { "M-D", 0x47 },
        { "D",   0x0C }, { "D",   0x4C }, 
        { "!D",  0x0D }, { "!D",  0x4D },
        { "D-1", 0x0E }, { "D-1", 0x4E },
        { "-D",  0x0F }, { "-D",  0x4F },
        { "D-A", 0x13 }, { "D-M", 0x53 },
        { "D|A", 0x15 }, { "D|M", 0x55 },
        { "D+1", 0x1F }, { "D+1", 0x5F }, 
        { "0",   0x2A }, { "0",   0x60 },
        { "A",   0x30 }, { "M",   0x70 },
        { "!A",  0x31 }, { "!M",  0x71 }, 
        { "A-1", 0x32 }, { "M-1", 0x72 },
        { "-A",  0x33 }, { "-M",  0x73 },
        { "A+1", 0x37 }, { "M+1", 0x77 },
        { "-1",  0x3A }, { "-1",  0x7A },
        { "1",   0x3F }, { "1",   0x7F }
};

char* instruction_bit16_to_asm(char* str, bit16_t* bit16)
{
    uint8_t is_c_instruction = bit16->bit[15];
    uint8_t is_a_instruction = ! is_c_instruction;

    /* A - instruction */
    if (is_a_instruction) {
        sprintf(str, "@%d # 0x%x", bit16_to_int16(bit16), bit16_to_int16(bit16));
        return str;
    }

    uint8_t jump_bits = 0;
            jump_bits |= bit16->bit[0] ? 0x1 : 0x0;
            jump_bits |= bit16->bit[1] ? 0x2 : 0x0;
            jump_bits |= bit16->bit[2] ? 0x4 : 0x0;
    char jump_bits_str[BIT16_SIZE];

    uint8_t dest_bits = 0;
            dest_bits |= bit16->bit[3] ? 0x1 : 0x0;
            dest_bits |= bit16->bit[4] ? 0x2 : 0x0;
            dest_bits |= bit16->bit[5] ? 0x4 : 0x0;
    char dest_bits_str[BIT16_SIZE];

    uint8_t comp_bits = 0;
            comp_bits |= bit16->bit[6]  ? 0x1 : 0x0;
            comp_bits |= bit16->bit[7]  ? 0x2 : 0x0;
            comp_bits |= bit16->bit[8]  ? 0x4 : 0x0;
            comp_bits |= bit16->bit[9]  ? 0x8 : 0x0;
            comp_bits |= bit16->bit[10] ? 0x10 : 0x0;
            comp_bits |= bit16->bit[11] ? 0x20 : 0x0;
            comp_bits |= bit16->bit[12] ? 0x40 : 0x0;
    char comp_bits_str[BIT16_SIZE];

    char* p_comp_bit_str = "???";
    int i;
    for (i=0; i < AR_SIZEOF(asm_a_comp_st); i++) {
	    if (comp_bits == asm_a_comp_st[i].bits) {
		p_comp_bit_str = asm_a_comp_st[i].name;
		break;
	    }
    }

    sprintf(str, "%s = %s ; %s # dest=%s, comp=%s, jump=%s",
        asm_a_dest_name[dest_bits], p_comp_bit_str, asm_a_jump_name[jump_bits],
        uint16_bits_to_string(dest_bits_str, dest_bits, 3),
        uint16_bits_to_string(comp_bits_str, comp_bits, 7),
        uint16_bits_to_string(jump_bits_str, jump_bits, 3));
    return str;
}

/*
 *         |<==============================*============================|
 *         |                               |       |f3                  |
 *         |==>|-\      +-----+            |    +--V--+   +----+        |
 * instr       |  |====>|reg_A|====>|      |===>|reg_D|==>|     \       |
 *  ==========>|-/      +--^--+     |           +-----+    \ ALU |      |    out_M
 *              ^          |f2      |                       |    |======|========> 
                |f1                 *============>|-\      /     |         write_M
   in_m                             |             |  |===>|     /        ^------->
 *  =============================================>|-/     +----+         |f6
 *                                  |              ^         ^
 *                                  |              |f4       |ccccc         addr_M
 *                                  *============================================>
 *                                  |               |f5
 *                                  |           +---V--+                     PC
 *                                  |==========>|cnt_PC|=========================>
 *                                              +------+
 */

	// A instruction:   0vvvvvvvvvvvvvvv        v => A
	// C instruction:   1xxaccccccdddjjj

/*      a==0 a==1   c c c c c c       dest  d d d    jump  j j j 
 *        0    0    1 0 1 0 1 0       null  0 0 0    null  0 0 0
 *        1    1    1 1 1 1 1 1        M    0 0 1    jgt   0 0 1 >0 (pos)
 *       -1   -1    1 1 1 0 1 0        D    0 1 0    jeq   0 1 0 =0
 *        D    D    0 0 1 1 0 0        A    1 0 0    jlt   1 0 0 <0 (neg) 
 *        A    M    1 1 0 0 0 0       DM    0 1 1    jge   0 1 1 >= 0  !neg
 *       !D   !D    0 0 1 1 0 1       AM    1 0 1    jne   1 0 1 !=0
 *       !A   !M    1 1 0 0 0 1       AD    1 1 0    jle   1 1 0 <=0 (!pos) 
 *       -D   -D    0 0 1 1 1 1       ADM   1 1 1    jmp   1 1 1 
 *       -A   -M    1 1 0 0 1 1
 *      D+1  D+1    0 1 1 1 1 1 
 *      A+1  M+1    1 1 0 1 1 1
 *      D-1  D-1    0 0 1 1 1 0
 *      A-1  M-1    1 1 0 0 1 0
 *      D+A  D+M    0 0 0 0 1 0 
 *      D-A  D-M    0 1 0 0 1 1
 *      A-D  M-D    0 0 0 1 1 1
 *      D&A  D&M    0 0 0 0 0 0 
 *      D|A  D|M    0 1 0 1 0 1
 */



void cpu_get_outputs_internal(cpu_t* cpu, bit16_t* alu_out, bit_t* alu_nz, bit_t* alu_ng, bit_t* wr_m, bit16_t* addr_m, bit16_t* pc)
{
    // A16, D16
    bit16_t a16 = reg16_get_output(&(cpu->reg_a16));
    bit16_t d16 = reg16_get_output(&(cpu->reg_d16));

    // PC16
    *pc = pc16_get_output(&(cpu->pc16));

    // Instruction = f(PC16)
    bit16_t instruction = rom32k_get_output(cpu->ext_rom, pc);

    // Memory addr
    *addr_m = a16;

    // Memory input = f(memory addr)
    bit16_t in_m = memory_get_output(cpu->ext_memory, addr_m);

    // ALU
	bit16_t alu_x, alu_y;
	/// x
	alu_x = d16;
	/// y
	bit_t sel_alu_y = instruction.bit[12];
	mux16(&a16, &in_m, sel_alu_y, &alu_y);
	/// control
    bit_t alu_zx = instruction.bit[11];
	bit_t alu_nx = instruction.bit[10];
	bit_t alu_zy = instruction.bit[9];
	bit_t alu_ny = instruction.bit[8];
	bit_t alu_f  = instruction.bit[7];
	bit_t alu_nf = instruction.bit[6];
	/// out
	alu(&alu_x, &alu_y, alu_zx, alu_nx, alu_zy, alu_ny, alu_f, alu_nf, alu_out, alu_nz, alu_ng);

    // Memory wr
    bit_t is_c_instruction = instruction.bit[15];
	gate_and(is_c_instruction, instruction.bit[3] /* d3: dest = M */, wr_m);
}

void cpu_get_outputs(cpu_t* cpu, bit16_t* out_m, bit_t* wr_m, bit16_t* addr_m, bit16_t* pc)
{
    bit16_t alu_out;
    bit_t alu_nz, alu_ng;
    cpu_get_outputs_internal(cpu, &alu_out, &alu_nz, &alu_ng, wr_m, addr_m, pc);
    *out_m = alu_out;
}

void cpu_run(cpu_t* cpu, bit_t reset)
{
    bit16_t   alu_out;
    bit_t alu_nz, alu_ng;
    bit_t      wr_m;
    bit16_t  addr_m;
    bit16_t    pc;

    // Get outputs
    cpu_get_outputs_internal(cpu, &alu_out, &alu_nz, &alu_ng, &wr_m /*isn't used*/, &addr_m /*isn't used*/, &pc);

    // Instruction = f(PC16)
    bit16_t instruction = rom32k_get_output(cpu->ext_rom, &pc);

    // Control = f(instruction)
	bit_t is_a_instruction;
	bit_t is_c_instruction = instruction.bit[15];
	gate_not(is_c_instruction, &is_a_instruction);

    // A16
    bit16_t a16 = reg16_get_output(&(cpu->reg_a16));
    bit_t sel_a = is_a_instruction;
    bit16_t in_a;
	mux16(&alu_out, &instruction, sel_a, &in_a);
    bit_t load_a;
    gate_or(is_a_instruction, instruction.bit[5]  /* d1: dest = A */, &load_a);

    // D16
	bit_t load_d;
    gate_and(is_c_instruction, instruction.bit[4] /* d2: dest = D */, &load_d);

    // PC16
    bit_t load_pc;
    bit_t jgt = instruction.bit[0];
	bit_t jeq = instruction.bit[1];
	bit_t jlt = instruction.bit[2];
	bit_t zero;
	gate_not(alu_nz, &zero);
	bit_t pos;
	gate_not(alu_ng, &pos);
	gate_and(pos, alu_nz, &pos);
	gate_and(jgt, pos,    &jgt);
	gate_and(jeq, zero,   &jeq);
	gate_and(jlt, alu_ng, &jlt);	
	or3_out1(jgt, jeq, jlt, &load_pc);
	gate_and(load_pc, is_c_instruction, &load_pc);

    /// RUN
    reg16_run(&(cpu->reg_a16), &in_a, load_a);
	reg16_run(&(cpu->reg_d16), &alu_out, load_d);
    pc16_run(&(cpu->pc16), reset, &a16, load_pc);
}

typedef struct {
	rom32k_t    rom32k;
	cpu_t       cpu;
	memory_t    memory;
} computer_t;

void computer_init(computer_t* computer)
{
	rom32k_init(&(computer->rom32k));
	memory_init(&(computer->memory));
	cpu_init(&(computer->cpu), &(computer->memory), &(computer->rom32k));	
}

void computer_run(computer_t* computer, bit_t reset)
{
    bit16_t pc;
    bit16_t out_m;
    bit16_t addr_m;
    bit_t   wr_m;

    // Get outputs
    cpu_get_outputs(&(computer->cpu), &out_m, &wr_m, &addr_m, &pc/*isn't used*/);

    // Run
    cpu_run(&(computer->cpu), reset);
    memory_run(&(computer->memory), &out_m, wr_m, &addr_m);
}

enum {
    ASM_OUT_FORMAT_STRIP = 1,
    ASM_OUT_FORMAT_SYMBOLS = 2,
    ASM_OUT_FORMAT_LITERALS = 3,
    ASM_OUT_FORMAT_BIN =4,
};

void asm_strip_line(char* line, char* strip_line)
{
    int i = 0;
    int j = 0;

    while (line[i] != '\0') {
        if (isspace(line[i])) {
            i++;
            continue;
        }
        if (line[i] == '/' && line[i+1] == '/') { // comment
            break;
        }
        strip_line[j++] = line[i++];
    }
    strip_line[j] = '\0';
}

enum {
    ASM_INSTRUCTION_FAIL             = 0,
    ASM_INSTRUCTION_C_DEST_COMP_JUMP = 1,
    ASM_INSTRUCTION_C_DEST_COMP      = 2,
    ASM_INSTRUCTION_C_COMP_JUMP      = 3,
    ASM_INSTRUCTION_C_COMP           = 4,
    ASM_INSTRUCTION_A_LITERAL        = 5,
    ASM_INSTRUCTION_A_SYMBOL         = 6,
    ASM_INSTRUCTION_L_SYMBOL         = 7
};

uint8_t asm_get_instruction_type(char* line)
{
    int n_at = 0; // @
    int n_lp = 0; // (
    int n_rp = 0; // )
    int n_e  = 0; // =
    int n_s  = 0; // ;
    
    int i = 0;
    while (line[i] != '\0') {
        if (line[i] == '@')
            n_at++;
        if (line[i] == '(')
            n_lp++;
        if (line[i] == ')')
            n_rp++;
        if (line[i] == '=')
            n_e++;
        if (line[i] == ';')
            n_s++;

        ++i;
    }
    printf("asm: @=%d, (=%d, )=%d, ==%d, :=%d\n", n_at, n_lp, n_rp, n_e, n_s);

    if (n_at == 1 && !n_lp && !n_rp && !n_e && !n_s) {
        if (isdigit(line[1]))
            return ASM_INSTRUCTION_A_LITERAL;
        return ASM_INSTRUCTION_A_SYMBOL;
    }

    if (!n_at && n_lp == 1 && n_rp == 1 && !n_e && !n_s) {
        if (isdigit(line[1]))
            return ASM_INSTRUCTION_FAIL; // only symbol
        return ASM_INSTRUCTION_L_SYMBOL; 
    }

    if (!n_at && !n_lp && !n_rp && n_e == 1 && n_s == 1)
        return ASM_INSTRUCTION_C_DEST_COMP_JUMP;
    if (!n_at && !n_lp && !n_rp && n_e == 1 && !n_s)
        return ASM_INSTRUCTION_C_DEST_COMP;
    if (!n_at && !n_lp && !n_rp && !n_e && n_s == 1)
        return ASM_INSTRUCTION_C_COMP_JUMP;
    if (!n_at && !n_lp && !n_rp && !n_e && !n_s)
        return ASM_INSTRUCTION_C_COMP;

    return ASM_INSTRUCTION_FAIL;
}

int asm_is_correct_symbol_name(char* symbol)
{
    if (strlen(symbol) == 0)
        return -1;
    int i = 0;
    uint8_t c;
    while((c=symbol[i++]) != '\0') {
        if (!isalpha(c) && !isdigit(c) && c != '_' && c != '.' && c != '$') {
           printf("asm: fail symbol char = %c\n", symbol[i-1]);
           return -1;
        }
    }
    return 0;
}

int asm_is_correct_literal(char* literal)
{
    if (strlen(literal) == 0) {
        printf("asm: empty literal\n");
        return -1;
    }
    int i = 0;
    uint8_t c;
    while((c=literal[i++]) != '\0') {
        if (!isdigit(c)) {
           printf("asm: fail literal char = %c\n", literal[i-1]);
           return -1;
        }
    }
    return 0;
}

int asm_get_symbol(char* strip_line, uint8_t instruction_type, char* symbol)
{
    int len = (int)strlen(strip_line);

    if (instruction_type == ASM_INSTRUCTION_A_SYMBOL) {
        if (strip_line[0] != '@' || len < 2) {
            printf("asm: fail A symbol l[0]=%c, len=%d\n", strip_line[0], len);
            return -1;
        }
        strcpy(symbol, &(strip_line[1]));
    }
    else if (instruction_type == ASM_INSTRUCTION_L_SYMBOL) {
        if (strip_line[0] != '(' || len < 3 || strip_line[len-1] != ')') {
            printf("asm: fail L symbol l[0]=%c, len=%d, l[-1]=%c\n", strip_line[0], len, strip_line[len-1]);
            return -1;
        }
        strcpy(symbol, &(strip_line[1]));
        symbol[len-1-1] = '\0'; // skip ')'
    }
    else {
        return -1;
    }

    if (asm_is_correct_symbol_name(symbol) < 0)
        return -1;

    return 0;
}

#define ASM_SYMBOL_MAX_LEN 256
#define ASM_SYMBOL_MAX_NUM 1024

typedef struct {
    char symbol_name[ASM_SYMBOL_MAX_LEN+1];
    int  symbol_value;

} asm_symbol_t;

typedef struct {
    asm_symbol_t symbols[ASM_SYMBOL_MAX_NUM];
    int n;
    int next_var_addr;
} asm_symbol_tbl_t;

void asm_add_name_symbol_tbl(asm_symbol_tbl_t* tbl, char* name)
{
    int i;
    for (i = 0; i < tbl->n; i++) {
        if (!strcmp(tbl->symbols[i].symbol_name, name))
            break;
    }
    // Add only if new name.
    if (i == tbl->n) {
        if (i < ASM_SYMBOL_MAX_NUM) {
            strcpy(tbl->symbols[i].symbol_name, name);
            tbl->symbols[i].symbol_value = -1;
        }
        ++(tbl->n);
    }
}

void asm_add_value_symbol_tbl(asm_symbol_tbl_t* tbl, char* name, int16_t value)
{
    int i;
    for (i = 0; i < tbl->n; i++) {
        if (!strcmp(tbl->symbols[i].symbol_name, name)) {
            if (tbl->symbols[i].symbol_value >= 0) {
                printf("asm: warning: rewrite symbol %s, old value = %d, new value %d\n", name, tbl->symbols[i].symbol_value, value);
            }
            tbl->symbols[i].symbol_value = value;
            return;
        }
    }
    // Add only if new name.
    if (i == tbl->n) {
        if (i < ASM_SYMBOL_MAX_NUM) {
            strcpy(tbl->symbols[i].symbol_name, name);
            tbl->symbols[i].symbol_value = value;
        }
        ++(tbl->n);
    }
}

int16_t asm_get_value_symbol_tbl(asm_symbol_tbl_t* tbl, char* name)
{
    int16_t value = -1;
    int i;
    for (i = 0; i < tbl->n; i++) {
        if (!strcmp(tbl->symbols[i].symbol_name, name)) {
            value = tbl->symbols[i].symbol_value;
            break;
        }
    }
    return value;
}

void asm_init_symbol_tbl(asm_symbol_tbl_t* tbl)
{
    tbl->n = 0;
    tbl->next_var_addr = 16;

    asm_add_value_symbol_tbl(tbl, "SP", 0);
    asm_add_value_symbol_tbl(tbl, "LCL", 1);
    asm_add_value_symbol_tbl(tbl, "ARG", 2);
    asm_add_value_symbol_tbl(tbl, "THIS", 3);
    asm_add_value_symbol_tbl(tbl, "THAT", 4);
    
    asm_add_value_symbol_tbl(tbl, "R0", 0);
    asm_add_value_symbol_tbl(tbl, "R1", 1);
    asm_add_value_symbol_tbl(tbl, "R2", 2);
    asm_add_value_symbol_tbl(tbl, "R3", 3);
    asm_add_value_symbol_tbl(tbl, "R4", 4);
    asm_add_value_symbol_tbl(tbl, "R5", 5);
    asm_add_value_symbol_tbl(tbl, "R6", 6);
    asm_add_value_symbol_tbl(tbl, "R7", 7);
    asm_add_value_symbol_tbl(tbl, "R8", 8);
    asm_add_value_symbol_tbl(tbl, "R9", 9);
    asm_add_value_symbol_tbl(tbl, "R10", 10);
    asm_add_value_symbol_tbl(tbl, "R11", 11);
    asm_add_value_symbol_tbl(tbl, "R12", 12);
    asm_add_value_symbol_tbl(tbl, "R13", 13);
    asm_add_value_symbol_tbl(tbl, "R14", 14);
    asm_add_value_symbol_tbl(tbl, "R15", 15);

    asm_add_value_symbol_tbl(tbl, "SCREEN", 16384);
    asm_add_value_symbol_tbl(tbl, "KBD", 24576);
}

void asm_set_variable_addr_symbol_tbl(asm_symbol_tbl_t* tbl)
{
    int i;
    for (i = 0; i < tbl->n; i++) {
        if (tbl->symbols[i].symbol_value < 0) {
            tbl->symbols[i].symbol_value = (tbl->next_var_addr)++;
        }
    }
}

#define ASM_LITERAL_MAX_LEN 256

int asm_decode(asm_symbol_tbl_t* symbol_tbl, char* strip_line, uint8_t instruction_type, uint16_t* code)
{
    int len = (int)strlen(strip_line);

    if (instruction_type == ASM_INSTRUCTION_A_LITERAL) {
        if (strip_line[0] != '@' || len < 2) {
            printf("asm: fail A literal l[0]=%c, len=%d\n", strip_line[0], len);
            return -1;
        }
        char literal[ASM_LITERAL_MAX_LEN+1];
        strcpy(literal, &(strip_line[1]));
        if (asm_is_correct_literal(literal) < 0)
            return -1;
        long literal_num = strtol(literal, NULL, 10);
        if (literal_num >= MEM_32K) {
            printf("asm: err literal is too much %ld\n", literal_num);
            return -1;
        }
        *code = (uint16_t)literal_num;
    }
    else if (instruction_type == ASM_INSTRUCTION_C_DEST_COMP_JUMP ||
             instruction_type == ASM_INSTRUCTION_C_DEST_COMP ||
             instruction_type == ASM_INSTRUCTION_C_COMP_JUMP ||
             instruction_type == ASM_INSTRUCTION_C_COMP) {
        *code = 0b1110000000000000;
        
        char dest[ASM_LITERAL_MAX_LEN+1] = "NULL";
        char comp[ASM_LITERAL_MAX_LEN+1] = "";
        char jump[ASM_LITERAL_MAX_LEN+1] = "NULL";
        if (instruction_type == ASM_INSTRUCTION_C_DEST_COMP_JUMP ||
            instruction_type == ASM_INSTRUCTION_C_DEST_COMP) {
             int i = 0;
             while (strip_line[i] != '\0') {
                if (strip_line[i] == '=')
                    break;
                i++;
             }
             strcpy(dest, strip_line);
             dest[i] = '\0';
        }
        if (instruction_type == ASM_INSTRUCTION_C_DEST_COMP_JUMP ||
            instruction_type == ASM_INSTRUCTION_C_COMP_JUMP) {
             int i = 0;
             while (strip_line[i] != '\0') {
                if (strip_line[i] == ';')
                    break;
                i++;
             }
             strcpy(jump, strip_line+i+1);
        }
        int i = 0;
        while (strip_line[i] != '\0') {
            if (strip_line[i] == '=')
                break;
            i++;
        }
        int j = 0;
        while (strip_line[j] != '\0') {
            if (strip_line[j] == ';')
                break;
            j++;
        }
        strcpy(comp, strip_line + ((i != len) ? (i+1) : 0));
        comp[j] = '\0';

        printf("asm: %s=%s;%s\n", dest, comp, jump);

        for (i = 0; i < AR_SIZEOF(asm_a_dest_st); i++)
            if (!strcmp(dest, asm_a_dest_st[i].name))
                break;
        if (i >= AR_SIZEOF(asm_a_dest_st)) {
            printf("asm: - err dest name %s\n", dest);
            return -1;
        }
        *code |= asm_a_dest_st[i].bits << 3;

        for (i = 0; i < AR_SIZEOF(asm_a_comp_st); i++)
            if (!strcmp(comp, asm_a_comp_st[i].name))
                break;
        if (i >= AR_SIZEOF(asm_a_comp_st)) {
            printf("asm: - err comp name %s\n", comp);
            return -1;
        }
        *code |= asm_a_comp_st[i].bits << 6;

        for (i = 0; i < AR_SIZEOF(asm_a_jump_name); i++)
            if (!strcmp(jump, asm_a_jump_name[i]))
                break;
        if (i >= AR_SIZEOF(asm_a_jump_name)) {
            printf("asm: - err jump name %s\n", jump);
            return -1;
        }
        *code |= ((uint16_t)i);
    }
    else {
        printf("asm: unknown intruction type\n");
        return -1;
    }

    return 0;
}

int asm_translate(char* asm_file_name, uint8_t out_format, asm_symbol_tbl_t* symbol_tbl, char* out_file_name)
{
    FILE* asm_file;
    FILE* out_file;

    asm_file = fopen(asm_file_name, "r");
    if (asm_file == NULL) {
        printf("asm: err to open asm file %s\n", asm_file_name);
        return -1;
    }

    out_file = fopen(out_file_name, "w");
    if (out_file == NULL) {
        printf("asm: err to create out file %s\n", out_file_name);
        fclose(asm_file);
        return -1;
    }

    #define ASM_LINE_MAX_LEN 1024
    char line[ASM_LINE_MAX_LEN+1];
    char strip_line[ASM_LINE_MAX_LEN+1];
    int n = 0;

    char symbol[ASM_SYMBOL_MAX_LEN+1];

    while (fgets(line, sizeof(line), asm_file) != NULL) {
        printf("asm: line %d <%s>, len = %d\n", n, line, (int)strlen(line));
        if (strlen(line) >= ASM_LINE_MAX_LEN) {
            printf("asm: -err line %d - str of line is too big = %d \n", n, (int)strlen(line));
            break;
        }
        asm_strip_line(line, strip_line);
        printf("asm: strip_line %d <%s>, len = %d\n", n, strip_line, (int)strlen(strip_line));
        if (strlen(strip_line) == 0)
            continue;
        if (out_format == ASM_OUT_FORMAT_STRIP) {
            fwrite(strip_line, 1, strlen(strip_line), out_file);
            fwrite("\n", 1, 1, out_file);
        }

        if (out_format == ASM_OUT_FORMAT_STRIP)
            continue;

        uint8_t instruction_type = asm_get_instruction_type(strip_line);
        printf("asm: line %d instruction type = %d\n", n, instruction_type);

        if (instruction_type == ASM_INSTRUCTION_FAIL) {
            printf("asm: -err line %d - fail instruction type\n", n);
            break;
        }

        if (instruction_type == ASM_INSTRUCTION_A_SYMBOL ||
            instruction_type == ASM_INSTRUCTION_L_SYMBOL) {

            if (asm_get_symbol(strip_line, instruction_type, symbol) < 0) {
                printf("asm: -err line %d - fail symbol\n", n);
                break;
            }

            int16_t symbol_addr = -1;
            if (instruction_type == ASM_INSTRUCTION_L_SYMBOL)
                symbol_addr = n;

            if (symbol_addr >= 0)
                asm_add_value_symbol_tbl(symbol_tbl, symbol, symbol_addr);
            else
                asm_add_name_symbol_tbl(symbol_tbl, symbol);

            if (out_format == ASM_OUT_FORMAT_SYMBOLS) {
                fwrite(strip_line, 1, strlen(strip_line), out_file);
                char asm_comment[ASM_LINE_MAX_LEN+1];
                sprintf(asm_comment, "    // ASM: line %d - symbol %s = %d\n", n, symbol, symbol_addr);
                fwrite(asm_comment, 1, strlen(asm_comment), out_file);
            }

            if (instruction_type == ASM_INSTRUCTION_L_SYMBOL)
                continue;

            if (out_format == ASM_OUT_FORMAT_LITERALS ||
                out_format == ASM_OUT_FORMAT_BIN) {
                symbol_addr = asm_get_value_symbol_tbl(symbol_tbl, symbol);
                if (symbol_addr < 0) {
                    printf("asm: -err line %d - unaddressed variable %s\n", n, symbol);
                    break;
                }
                if (out_format == ASM_OUT_FORMAT_LITERALS) {
                    char asm_literal[ASM_LINE_MAX_LEN+1];
                    sprintf(asm_literal, "@%d\n", symbol_addr);
                    fwrite(asm_literal, 1, strlen(asm_literal), out_file);
                }
                char binary_str[16+1+1];
                uint16_t u = symbol_addr;
                int i;
                for (i=0; i<16; i++) {
                    binary_str[16-1-i] =  (u & 0x01) ? '1' : '0';
                    u >>= 1;
                }
                binary_str[16] = '\n';
                binary_str[17] = '\0';

                if (out_format == ASM_OUT_FORMAT_BIN) {
                    fwrite(binary_str, 1, strlen(binary_str), out_file);
                }
            }
        }
        else {
            if (out_format == ASM_OUT_FORMAT_SYMBOLS ||
                out_format == ASM_OUT_FORMAT_LITERALS) 
            {
                fwrite(strip_line, 1, strlen(strip_line), out_file);
                fwrite("\n", 1, 1, out_file);
            }
        }

        n++;
        
        if (out_format == ASM_OUT_FORMAT_SYMBOLS || out_format == ASM_OUT_FORMAT_LITERALS)
            continue;

        if (out_format == ASM_OUT_FORMAT_BIN && instruction_type == ASM_INSTRUCTION_A_SYMBOL)
            continue;

        uint16_t code;
        if (asm_decode(symbol_tbl, strip_line, instruction_type, &code) < 0)
            break;
        char binary_str[16+1+1];
        uint16_t u = code;
        int i;
        for (i=0; i<16; i++) {
            binary_str[16-1-i] =  (u & 0x01) ? '1' : '0';
            u >>= 1;
        }
        binary_str[16] = '\n';
        binary_str[17] = '\0';

        if (out_format == ASM_OUT_FORMAT_BIN) {
            fwrite(binary_str, 1, strlen(binary_str), out_file);
        }
    }
    
    if (out_format == ASM_OUT_FORMAT_SYMBOLS) {
        asm_set_variable_addr_symbol_tbl(symbol_tbl);

        int i;
        for (i = 0; i < symbol_tbl->n; i++) {
            char asm_comment[ASM_LINE_MAX_LEN+1] = {'\0'};
            sprintf(asm_comment, "// ASM: stbl %d. %s = %d (0x%.4x)\n", i,
                    symbol_tbl->symbols[i].symbol_name, symbol_tbl->symbols[i].symbol_value,
                    symbol_tbl->symbols[i].symbol_value);
            fwrite(asm_comment, 1, strlen(asm_comment), out_file);
        }
    }

    fclose(asm_file);
    fclose(out_file);
    return 0;
}

enum {
    VM_OUT_FORMAT_STRIP = 1,
    VM_OUT_FORMAT_ASM = 2
};

void vm_strip_line(char* line, char* strip_line)
{
    int i = 0;
    int j = 0;
    uint8_t is_space = 0;

    while (line[i] != '\0') {
        if (line[i] == '/' && line[i+1] == '/') { // comment
            break;
        }
        if (isspace(line[i])) {
            i++;
            is_space = 1;
            continue;
        }
        else {
            if (is_space) {
                if (j)
                    strip_line[j++] = ' ';
                is_space = 0;
            }
        }

        strip_line[j++] = line[i++];
    }
    strip_line[j] = '\0';
}

#define VM_COMMAND_LINE_MAX 256

char* vm_get_command_line(char* line, char* command)
{
    int i = 0;
    int j = 0;
    uint8_t c;
    while ((c = line[i++]) != '\0') {
        if (!isalpha(c) && !isdigit(c) && c != '-') {
            if (c != ' ') {
                printf("vm: - error cmd has symbol %c\n", c);
                return NULL;
            }
            if (j == 0) {
                printf("vm: - error emtpy cmd\n");
                return NULL;
            }
            break;
        }
        if (j >= VM_COMMAND_LINE_MAX) {
            printf("vm: -error very long command >= %d\n", VM_COMMAND_LINE_MAX);
            return NULL;
        }
        command[j++] = c;
    }
    command[j] = '\0';
    return (line + i);
}

enum {
    VM_INSTRUCTION_FAIL             = 0,
    VM_INSTRUCTION_FUNCTION         = 1,
    VM_INSTRUCTION_PUSH             = 2,
    VM_INSTRUCTION_LT               = 3,
    VM_INSTRUCTION_IF_GOTO          = 4,
    VM_INSTRUCTION_GOTO             = 5,
    VM_INSTRUCTION_LABEL            = 6,
    VM_INSTRUCTION_POP              = 7,
    VM_INSTRUCTION_RETURN           = 8,
    VM_INSTRUCTION_SUB              = 9,
    VM_INSTRUCTION_CALL             = 10,
    VM_INSTRUCTION_ADD              = 11,
};

struct {
    char* name;
    uint8_t type;
} vm_command_st[] = {
    {"function",     VM_INSTRUCTION_FUNCTION},
    {"push",         VM_INSTRUCTION_PUSH},
    {"lt",           VM_INSTRUCTION_LT},
    {"if-goto",      VM_INSTRUCTION_IF_GOTO},
    {"goto",         VM_INSTRUCTION_GOTO},
    {"label",        VM_INSTRUCTION_LABEL},
    {"pop",          VM_INSTRUCTION_POP},
    {"return",       VM_INSTRUCTION_RETURN},
    {"sub",          VM_INSTRUCTION_SUB},
    {"call",         VM_INSTRUCTION_CALL},
    {"add",          VM_INSTRUCTION_ADD},
};

uint8_t vm_get_instruction_type(char* line, char** next)
{
    char command_line[VM_COMMAND_LINE_MAX+1];

    *next = vm_get_command_line(line, command_line);
    if (*next == NULL)
        return VM_INSTRUCTION_FAIL;

    printf("vm: command_line = %s; next = <%s>\n", command_line, *next);

    uint8_t command_type = VM_INSTRUCTION_FAIL;
    int i;
    for(i = 0; i < AR_SIZEOF(vm_command_st); i++) {
        if (!strcmp(command_line, vm_command_st[i].name)) {
            command_type = vm_command_st[i].type;
            break;
        }
    }

    return command_type;
}

char* vm_get_label_line(char* line, char* label)
{
    int i = 0;
    int j = 0;
    uint8_t c;
    while ((c=line[i++]) != '\0') {
        if (isspace(c)) {
            if (j == 0) {
                printf("vm: - error empty label\n");
                return NULL;
            }
            break;
        }
        label[j++] = c;
    }
    label[j] = '\0';

    if (asm_is_correct_symbol_name(label) < 0) {
        return NULL;
    }

    return (line + i);
}

char* vm_get_literal(char* line, char* literal)
{
    int i = 0;
    int j = 0;
    uint8_t c;
    while ((c=line[i++]) != '\0') {
        if (!isdigit(c)) {
            if (c != ' ') {
                printf("vm; - error character %c in literal\n", c);
                return NULL;
            }
            if (j == 0) {
                printf("vm: - error empty literal\n");
                return NULL;
            }
            break;
        }
        literal[j++] = c;
    }
    literal[j] = '\0';

    if (asm_is_correct_literal(literal) < 0) {
        return NULL;
    }

    return (line + i);
}

int vm_write_push_d_and_inc_sp(char* asm_lines)
{
    sprintf(asm_lines, 
                // *(*SP) = D; SP++
                /// *(*SP) = D
                "@SP\n"
                "A=M\n"
                "M=D\n"
                /// SP++
                "@SP\n"
                "M=M+1\n");
    return 0;
}

int vm_write_pop_d_and_dec_sp(char* asm_lines)
{
    sprintf(asm_lines, 
                /// D = *(SP-1); SP = SP-1
                "@SP\n"
                "AM=M-1\n"
                "D=M\n");
    return 0;
}
    

int vm_write_push_reg_index(char* reg, uint16_t n, char* asm_lines)
{
    // *(*SP) = *(*(reg) + n)
    /// D = *(*(reg) + n)
    if (n > 1) {
        sprintf(asm_lines, 
                    "@%u\n"
                    "D=A\n", n);
        asm_lines += strlen(asm_lines);
    }
    sprintf(asm_lines, 
                    "@%s\n", reg);
    asm_lines += strlen(asm_lines);
    if (n == 0) {
        sprintf(asm_lines, 
                    "A=M\n");
    }
    else if (n == 1) {
        sprintf(asm_lines, 
                    "A=M+1\n");
    }
    else {
        sprintf(asm_lines, 
                    "A=D+M\n");
    }
    asm_lines += strlen(asm_lines);
    sprintf(asm_lines, 
                    "D=M\n");
    vm_write_push_d_and_inc_sp(asm_lines + strlen(asm_lines));
    return 0;
}

int vm_write_pop_reg_index(char* reg, uint16_t n, char* asm_lines)
{
    // *(*(reg) + n) = *(*(SP - 1)); SP--
    ////  R13 = *(reg) + n
    if (n > 1) {
        sprintf(asm_lines, 
                    "@%u\n"
                    "D=A\n", n);
        asm_lines += strlen(asm_lines);
    }
    sprintf(asm_lines, 
                    "@%s\n", reg);
    asm_lines += strlen(asm_lines);
    if (n == 0) {
        sprintf(asm_lines, 
                    "D=M\n");
    }
    else if (n == 1) {
        sprintf(asm_lines, 
                    "D=M+1\n");
    }
    else {
        sprintf(asm_lines, 
                    "D=D+M\n");
    }
    asm_lines += strlen(asm_lines);
    sprintf(asm_lines, 
                    "@R13\n"
                    "M=D\n");
    ////  D = *(SP-1); SP = SP-1
    vm_write_pop_d_and_dec_sp(asm_lines + strlen(asm_lines));
    asm_lines += strlen(asm_lines);
    ////  *(R13) = D
    sprintf(asm_lines, 
                    "@R13\n"
                    "A=M\n"
                    "M=D\n");
    return 0;
}

int vm_write_push(char* segment, char* literal, char* asm_lines)
{
    uint16_t literal_num = (uint16_t)strtol(literal, NULL, 10);

    if (!strcmp(segment, "constant")) {
        // D = a
        sprintf(asm_lines, 
                    // D = a
                    "@%u\n"
                    "D=A\n",
                    literal_num);
        vm_write_push_d_and_inc_sp(asm_lines + strlen(asm_lines));
        return 0;
    }
    else if (!strcmp(segment, "argument")) {
        if (vm_write_push_reg_index("ARG", literal_num, asm_lines) < 0)
            return -1;
        return 0;
    }
    
    printf("vm: - error push unknow segment %s\n", segment);
    return -1; 
}

int vm_write_pop(char* segment, char* literal, char* asm_lines)
{
    uint16_t literal_num = (uint16_t)strtol(literal, NULL, 10);

    if (!strcmp(segment, "argument")) {
        if (vm_write_pop_reg_index("ARG", literal_num, asm_lines) < 0)
            return -1;
        return 0;
    }
    
    printf("vm: - error pop unknow segment %s\n", segment);
    return -1; 
}

int vm_write_function(char* label, char* literal, char* asm_lines)
{
    int16_t n = (int16_t)strtol(literal, NULL, 10);

    sprintf(asm_lines, "(%s)\n", label);
    int i;
    for (i = 0; i < n; i++) {
        if (vm_write_push("constant", "0", (asm_lines + strlen(asm_lines))) < 0)
            return -1;
    }
    return 0;
}

int vm_write_push_reg(char* reg, char* asm_lines)
{
    sprintf(asm_lines,
                "@%s\n"
                "D=M\n", reg);
    vm_write_push_d_and_inc_sp(asm_lines + strlen(asm_lines));
    return 0;
}

uint16_t s_returtn_number = 0;

int vm_write_call(char* label, char* literal, char* asm_lines)
{
    int16_t n = (int16_t)strtol(literal, NULL, 10);

    //push returnAddress	- (generate a label and pushes it to the stack)
    sprintf(asm_lines,
                "@%s$ret.%u\n"
                "D=A\n", label, s_returtn_number);
    vm_write_push_d_and_inc_sp(asm_lines + strlen(asm_lines));
    
    vm_write_push_reg("LCL", asm_lines + strlen(asm_lines));
    vm_write_push_reg("ARG", asm_lines + strlen(asm_lines));
    vm_write_push_reg("THIS", asm_lines + strlen(asm_lines));
    vm_write_push_reg("THAT", asm_lines + strlen(asm_lines));

    // ARG = SP - 5 - nArgs
    sprintf(asm_lines + strlen(asm_lines),
                "@%d\n"
                "D=A\n"
                "@SP\n"
                "D=M-D\n"
                "@ARG\n"
                "M=D\n"
                "@SP\n"
                "D=M\n"
                // LCL = SP
                "@LCL\n"
                "M=D\n",
                n + 5);
                
    // goto f
    sprintf(asm_lines + strlen(asm_lines),
                "@%s\n"
                "0;JMP\n", label);

    // (return address)
    sprintf(asm_lines + strlen(asm_lines),
                "(%s$ret.%u)\n", label, s_returtn_number);

    s_returtn_number++;

    return 0;
}

uint16_t s_label_num = 0;

int vm_write_compare(char* cmd, char* asm_lines)
{
    vm_write_pop_d_and_dec_sp(asm_lines);
    sprintf(asm_lines + strlen(asm_lines), 
                "A=A-1\n"
                "D=M-D\n"
                "@%s_TRUE_%u\n"
                "D;J%s\n"
                "@SP\n"
                "A=M-1\n"
                "M=0\n"
                "@%s_END_%u\n"
                "0;JMP\n"
                "(%s_TRUE_%u)\n"
                "@SP\n"
                "A=M-1\n"
                "M=-1\n"
                "(%s_END_%u)\n",
                cmd, s_label_num,
                cmd,
                cmd, s_label_num,
                cmd, s_label_num,
                cmd, s_label_num);
    s_label_num++;
    return 0;
}

int vm_write_goto(char* func_name, char* label, char* asm_lines)
{
    sprintf(asm_lines,
                "@%s$%s\n"
                "0;JMP\n",
                func_name, label);
    return 0;
}

int vm_write_if_goto(char* func_name, char* label, char* asm_lines)
{
    vm_write_pop_d_and_dec_sp(asm_lines);
    sprintf(asm_lines + strlen(asm_lines),
                "@%s$%s\n"
                "D;JNE\n",
                func_name, label);
    return 0;
}

int vm_write_binary_operator(char* operation, char* asm_lines)
{
    char comp[] = "X*Y";

    sprintf(comp, (!strcmp(operation, "-")) ? "M%sD" : "D%sM", operation);

    /// D = *(SP-1); SP = SP-1
    vm_write_pop_d_and_dec_sp(asm_lines);
    sprintf(asm_lines + strlen(asm_lines),
                "A=A-1\n"
                "M=%s\n",
                comp);

    return 0;
}

int vm_write_label(char* func_name, char* label, char* asm_lines)
{
    sprintf(asm_lines,
                "(%s$%s)\n",
                func_name, label);
    return 0;
}

int vm_write_return(char* asm_lines)
{
    sprintf(asm_lines,
                // local frame = LCL
                "@LCL\n"
                "D=M\n"
                "@R14\n"
                "M=D\n"
                // local retAddr = *(frame-5)
                "@5\n"
                "D=A\n"
                "@R14\n"
                "A=M-D\n"
                "D=M\n"
                "@R15\n"
                "M=D\n");
    //pop argument 0
    if (vm_write_pop_reg_index("ARG", 0, asm_lines + strlen(asm_lines)) < 0)
            return -1;
    asm_lines += strlen(asm_lines);
    sprintf(asm_lines,
                // SP = ARG+1
                "@ARG\n"
                "D=M+1\n"
                "@SP\n"
                "M=D\n"
                // THAT = *(--frame))
                "@R14\n"
                "AM=M-1\n"
                "D=M\n"
                "@THAT\n"
                "M=D\n"
                // THIS = *(--frame)
                "@R14\n"
                "AM=M-1\n"
                "D=M\n"
                "@THIS\n"
                "M=D\n"
                // ARG = *(--frame))
                "@R14\n"
                "AM=M-1\n"
                "D=M\n"
                "@ARG\n"
                "M=D\n"
                // LCL = *(--frame))
                "@R14\n"
                "AM=M-1\n"
                "D=M\n"
                "@LCL\n"
                "M=D\n"
                // goto retAddr
                "@R15\n"
                "A=M\n"
                "0;JMP\n");
    return 0;
}

int vm_write_asm(uint8_t instruction_type, char** next, char* func_name, char* asm_lines)
{
    char label[ASM_SYMBOL_MAX_LEN+1];
    char literal[ASM_LITERAL_MAX_LEN+1];

    if (instruction_type == VM_INSTRUCTION_FUNCTION) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: function name = %s, next = <%s>\n", label, *next);

        strcpy(func_name, label);
        
        *next = vm_get_literal(*next, literal);
        if (*next == NULL)
            return -1;
        printf("vm: function n param = %s, next = <%s>\n", literal, *next);

        if (vm_write_function(label, literal, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_PUSH) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: push segment = %s, next = <%s>\n", label, *next);
        
        *next = vm_get_literal(*next, literal);
        if (*next == NULL)
            return -1;
        printf("vm: push n = %s, next = <%s>\n", literal, *next);

        if (vm_write_push(label, literal, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_LT) {
        vm_write_compare("LT", asm_lines);
    }
    if (instruction_type == VM_INSTRUCTION_IF_GOTO) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: if-goto label = %s, next = <%s>\n", label, *next);

        if (vm_write_if_goto(func_name, label, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_GOTO) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: goto label = %s, next = <%s>\n", label, *next);

        if (vm_write_goto(func_name, label, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_LABEL) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: label = %s, next = <%s>\n", label, *next);

        if (vm_write_label(func_name, label, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_POP) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: pop segment = %s, next = <%s>\n", label, *next);
        
        *next = vm_get_literal(*next, literal);
        if (*next == NULL)
            return -1;
        printf("vm: pop n = %s, next = <%s>\n", literal, *next);

        if (vm_write_pop(label, literal, asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_RETURN) {
        if (vm_write_return(asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_SUB) {
        if (vm_write_binary_operator("-", asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_ADD) {
        if (vm_write_binary_operator("+", asm_lines) < 0)
            return -1;
    }
    if (instruction_type == VM_INSTRUCTION_CALL) {
        *next = vm_get_label_line(*next, label);
        if (*next == NULL)
            return -1;
        printf("vm: call function name = %s, next = <%s>\n", label, *next);

        *next = vm_get_literal(*next, literal);
        if (*next == NULL)
            return -1;
        printf("vm: call function n args = %s, next = <%s>\n", literal, *next);

        if (vm_write_call(label, literal, asm_lines) < 0)
            return -1;
    }

    return 0;
}

int vm_translate(char* vm_file_name, uint8_t out_format, char* out_file_name)
{
    FILE* vm_file;
    FILE* out_file;

    vm_file = fopen(vm_file_name, "r");
    if (vm_file == NULL) {
        printf("vm: err to open asm file %s\n", vm_file_name);
        return -1;
    }

    out_file = fopen(out_file_name, "w");
    if (out_file == NULL) {
        printf("vm: err to create out file %s\n", out_file_name);
        fclose(vm_file);
        return -1;
    }

    #define VM_LINE_MAX_LEN 1024
    char line[VM_LINE_MAX_LEN+1];
    char strip_line[VM_LINE_MAX_LEN+1];
    int n = 0;

    #define VM_ASM_LINES_MAX 2048
    char asm_lines[VM_ASM_LINES_MAX+1];

    char func_name[ASM_SYMBOL_MAX_LEN+1] = "";

    if (out_format == VM_OUT_FORMAT_ASM) {
          sprintf(asm_lines, "// asm: bootstrap\n"
                             "/// SP = 256\n"
                             "@256\n"
                             "D=A\n"
                             "@SP\n"
                             "M=D\n"
                             "/// call function Sys.init 0\n");
          vm_write_call("Sys.init", "0", asm_lines + strlen(asm_lines));
          fwrite(asm_lines, 1, strlen(asm_lines), out_file);
    }

    while (fgets(line, sizeof(line), vm_file) != NULL) {
        printf("vm: line %d <%s>, len = %d\n", n, line, (int)strlen(line));
        if (strlen(line) >= VM_LINE_MAX_LEN) {
            printf("vm: -err line %d - str of line is too big = %d \n", n, (int)strlen(line));
            break;
        }
        vm_strip_line(line, strip_line);
        printf("vm: strip_line %d <%s>, len = %d\n", n, strip_line, (int)strlen(strip_line));
        if (strlen(strip_line) == 0)
            continue;
        if (out_format == VM_OUT_FORMAT_STRIP) {
            fwrite(strip_line, 1, strlen(strip_line), out_file);
            fwrite("\n", 1, 1, out_file);
        }

        if (out_format == VM_OUT_FORMAT_STRIP)
            continue;

        char* next;
        uint8_t instruction_type = vm_get_instruction_type(strip_line, &next);
        printf("vm: line %d instruction type = %d\n", n, instruction_type);


        if (instruction_type == VM_INSTRUCTION_FAIL) {
            printf("vm: -err line %d - fail instruction type\n", n);
            break;
        }

        if (out_format == VM_OUT_FORMAT_ASM) {
            char str_comment[VM_LINE_MAX_LEN+1+256];
            sprintf(str_comment, "// vm: line %d <%s>\n($vm.line.%u)\n", n, strip_line, n);
            fwrite(str_comment, 1, strlen(str_comment), out_file);
        }

        n++;

        if (vm_write_asm(instruction_type, &next, func_name, asm_lines) < 0) {
            break;
        }
        if (out_format == VM_OUT_FORMAT_ASM) {
          fwrite(asm_lines, 1, strlen(asm_lines), out_file);
        }
    }

    fclose(vm_file);
    fclose(out_file);
    return 0;
}

enum {
    COMPILER_OUT_FORMAT_STRIP            = 1,
    COMPILER_OUT_FORMAT_XML_TOKENS       = 2,
    COMPILER_OUT_FORMAT_XML_PARSER       = 3,
    COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE = 4,
    COMPILER_OUT_FORMAT_VM               = 5,
};

#define COMPILER_LINE_MAX_LEN 1024

int compiler_strip_line(char* line, uint8_t* is_in_comment, char* strip_line)
{
    char strip_vm_line[COMPILER_LINE_MAX_LEN+1];
    vm_strip_line(line, strip_vm_line);
    printf("compiler: vm_strip_line <%s> len = %u\n", strip_vm_line, (unsigned)strlen(strip_vm_line));


    int i = 0;
    int j = 0;
    uint8_t c;

    while ((c = strip_vm_line[i]) != '\0') {
        if (*is_in_comment) {
            if (c == '*' && strip_vm_line[i+1] == '/') {
                if (!(*is_in_comment))
                    return -1;
                *is_in_comment = 0;
                i++;
            }
            i++;
            continue;
        }

        if (c == '/' && strip_vm_line[i+1] == '*') {
            *is_in_comment = 1;
            i+=2;
            continue;
        }

        strip_line[j++] = strip_vm_line[i++];
    }
    strip_line[j] = '\0';
    return 0;
}

uint8_t compiler_symbols[] = {
    '{', '}', '(', ')', '[', ']', '.', ',', ';',
    '+', '-', '*', '/', '&', '|', '<', '>', '=', '~'
};

int compiler_is_symbol(uint8_t c) {
    int i;
    for (i = 0; i < AR_SIZEOF(compiler_symbols); i++) {
        if (c == compiler_symbols[i])
            return 1;
    }
    return 0;
}

int compiler_get_token(char** line, char* token)
{
    uint8_t is_in_str_literal = 0;
    int i = 0;
    int j = 0;
    uint8_t c;
    while((c=(*line)[i]) != '\0') {
        if (c == '\"') {
            if (is_in_str_literal) {
                token[j++];
                i++;
                break;
            }
            is_in_str_literal = 1;
        }
        else if (is_in_str_literal) {
            // nothing
        }
        else if (c == ' ') {
            if (j == 0) {
                printf("compiler: empty token\n");
                return -1;
            }
            i++;
            break;
        }
        else if (compiler_is_symbol(c)) {
            if (j == 0) {
                token[j++] = c;
                i++;
                if ((*line)[i] == ' ')
                    i++;
            }
            break;
        }
        else if (!isalpha(c) && !isdigit(c) && c != '_') {
            printf("compiler: unknow symbol %c\n", c);
            return -1;
        }
        token[j++] = c;
        i++;
    }
    token[j] = '\0';
    *line += i;
    return 0;
}

int compiler_is_string_const(char* token)
{
    if (strlen(token) == 0) 
        return 0;

    return (token[0] == '\"') ? 1 : 0;
}

void compiler_skip_string_quotes(char* token)
{
    int len = strlen(token);
    if (len <= 2) {
        token[0] = '\0';
    }
        
    int i;
    for (i = 1; i < len; i++)
        token[i-1] = token[i];
    token[len-2] = '\0';
}

int compiler_is_int_const(char* token)
{
    if (strlen(token) == 0) 
        return 0;

    int i = 0;
    uint8_t c;
    while((c=token[i++]) != '\0') {
        if (!isdigit(c))
           return 0;
    }
    return 1;
}

int compiler_is_identifier(char* token)
{
    if (strlen(token) == 0) 
        return 0;

    if (isdigit(token[0]))
        return 0;

    int i = 0;
    uint8_t c;
    while((c=token[i++]) != '\0') {
        if (!isdigit(c) && !isalpha(c) && (c != '_'))
           return 0;
    }
    return 1;
}

enum {
    COMPILER_KEYWORD_CLASS           = 0,
    COMPILER_KEYWORD_CONSTRUCTOR     = 1,
    COMPILER_KEYWORD_FUNCTION        = 2,
    COMPILER_KEYWORD_METHOD          = 3,
    COMPILER_KEYWORD_FIELD           = 4,
    COMPILER_KEYWORD_STATIC          = 5,
    COMPILER_KEYWORD_VAR             = 6,
    COMPILER_KEYWORD_INT             = 7,
    COMPILER_KEYWORD_CHAR            = 8,
    COMPILER_KEYWORD_BOOLEAN         = 9,
    COMPILER_KEYWORD_VOID            = 10,
    COMPILER_KEYWORD_TRUE            = 11,
    COMPILER_KEYWORD_FALSE           = 12,
    COMPILER_KEYWORD_NULL            = 13,
    COMPILER_KEYWORD_THIS            = 14,
    COMPILER_KEYWORD_LET             = 15,
    COMPILER_KEYWORD_DO              = 16,
    COMPILER_KEYWORD_IF              = 17,
    COMPILER_KEYWORD_ELSE            = 18,
    COMPILER_KEYWORD_WHILE           = 19,
    COMPILER_KEYWORD_RETURN          = 20,
};

char* compiler_keywords[] = {
    "class",
    "constructor",
    "function",
    "method",
    "field",
    "static",
    "var",
    "int",
    "char",
    "boolean",
    "void",
    "true",
    "false",
    "null",
    "this",
    "let",
    "do",
    "if",
    "else",
    "while",
    "return",
};

int compiler_is_keyword(char* token)
{
    int i = 0;
    for (i = 0; i < AR_SIZEOF(compiler_keywords); i++) {
        if (!strcmp(token, compiler_keywords[i]))
            return 1;
    }
    return 0;
}

uint8_t compiler_get_keyword_type(char* token) {
    int i = 0;
    for (i = 0; i < AR_SIZEOF(compiler_keywords); i++) {
        if (!strcmp(token, compiler_keywords[i]))
            return (uint8_t)i;
    }
    return (uint8_t)i;
}

enum {
    COMPILER_TOKEN_NAME_UNKNOWN      = 0,
    COMPILER_TOKEN_NAME_KEYWORD      = 1,
    COMPILER_TOKEN_NAME_SYMBOL       = 2,
    COMPILER_TOKEN_NAME_IDENTIFIER   = 3,
    COMPILER_TOKEN_NAME_INT_CONST    = 4,
    COMPILER_TOKEN_NAME_STRING_CONST = 5,
    COMPILER_TOKEN_NAME_END
};
char* compiler_token_name[] = {
    "unknown",   // COMPILER_TOKEN_NAME_UNKNOWN
    "keyword",
    "symbol",
    "identifier",
    "intConst",
    "stringConst",
};

#define COMPILER_TOKEN_MAX_NUM 1024
#define COMPILER_TOKEN_MAX_LEN 1024
typedef struct {
    uint8_t     type;
    char        token[COMPILER_TOKEN_MAX_LEN+1];
} token_t;

typedef struct {
    token_t tokens[COMPILER_TOKEN_MAX_NUM];
    int     num;
} tokens_t;

void compiler_init_tokens(tokens_t* tokens) {
    tokens->num = 0;
}

int compiler_add_to_tokens(tokens_t* tokens, uint8_t type, char* token) {
    int i = tokens->num;

    if (i >= COMPILER_TOKEN_MAX_NUM)
        return -1;

    tokens->tokens[i].type = type;
    strcpy(tokens->tokens[i].token, token);

    ++(tokens->num);
    return 0;
}

token_t* compiler_get_token(tokens_t* tokens, int n) {
    if (n >= tokens->num)
        return NULL;

    return &(tokens->tokens[n]);
}

uint8_t compiler_get_token_type(char* token) {
    uint8_t token_type = COMPILER_TOKEN_NAME_UNKNOWN;

    if (strlen(token) == 1 && compiler_is_symbol(token[0]))
        token_type = COMPILER_TOKEN_NAME_SYMBOL;
    else if (compiler_is_int_const(token))
        token_type = COMPILER_TOKEN_NAME_INT_CONST;
    else if (compiler_is_string_const(token))
        token_type = COMPILER_TOKEN_NAME_STRING_CONST;
    else if (compiler_is_keyword(token))
        token_type = COMPILER_TOKEN_NAME_KEYWORD;
    else if (compiler_is_identifier(token))
        token_type = COMPILER_TOKEN_NAME_IDENTIFIER;
    
    return token_type;
}

uint8_t compiler_is_token_keyword_const(token_t* token)
{
    if (COMPILER_TOKEN_NAME_KEYWORD != token->type)
        return 0;

    uint8_t type = compiler_get_keyword_type(token->token);
    switch (type) {
        case COMPILER_KEYWORD_TRUE:
        case COMPILER_KEYWORD_FALSE:
        case COMPILER_KEYWORD_NULL:
        case COMPILER_KEYWORD_THIS:
            return 1;
        
        default:
            break;
    }

    return 0;
}

uint8_t compiler_is_token_op(token_t* token)
{
    if (COMPILER_TOKEN_NAME_SYMBOL != token->type)
        return 0;

    if ((!strcmp("+", token->token)) ||
        (!strcmp("-", token->token)) ||
        (!strcmp("*", token->token)) ||
        (!strcmp("/", token->token)) ||
        (!strcmp("&", token->token)) ||
        (!strcmp("|", token->token)) ||
        (!strcmp("<", token->token)) ||
        (!strcmp(">", token->token)) ||
        (!strcmp("=", token->token))) {
        return 1;
    }

    return 0;
}

void compiler_xml_print_start_tag(char* tag, int n_tabs, FILE* out_file) {

    char token_line[1024+1];
    int i;
    //
    token_line[0] = '\0';
    for (i = 0; i < n_tabs; i++)
        strcat(token_line, "    ");
    strcat(token_line, "<");
    strcat(token_line, tag);
    strcat(token_line, ">\n");
    fprintf(out_file, "%s", token_line);
}

void compiler_xml_print_end_tag(char* tag, int n_tabs, FILE* out_file) {

    char token_line[1024+1];
    int i;
    //
    token_line[0] = '\0';
    for (i = 0; i < n_tabs; i++)
        strcat(token_line, "    ");
    strcat(token_line, "</");
    strcat(token_line, tag);
    strcat(token_line, ">\n");
    fprintf(out_file, "%s", token_line);
}

void compiler_xml_print_token(token_t* token, int n_tabs, FILE* out_file)
{
    char token_line[1024+1];
    int i;
    //
    token_line[0] = '\0';
    for (i = 0; i < n_tabs; i++)
        strcat(token_line, "    ");
    strcat(token_line, "<");
    strcat(token_line, compiler_token_name[token->type]);
    strcat(token_line, "> ");
    strcat(token_line, (!strcmp(token->token, "<")) ? "&lt;" : token->token);
    strcat(token_line, " </");
    strcat(token_line, compiler_token_name[token->type]);
    strcat(token_line, ">\n");
    //
    fprintf(out_file, "%s", token_line);
}


enum compiler_parser_elem_type {
    COMPILER_PARSER_ELEM_CLASS                  = 0,
    COMPILER_PARSER_ELEM_CLASS_VAR_DEC          = 1,
    COMPILER_PARSER_ELEM_PARAMETER_LIST         = 2,
    COMPILER_PARSER_ELEM_VAR_DEC                = 3,
    COMPILER_PARSER_ELEM_EXPRESSION_LIST        = 4,
    COMPILER_PARSER_ELEM_RETURN_STATEMENT       = 5,
    COMPILER_PARSER_ELEM_DO_STATEMENT           = 6,
    COMPILER_PARSER_ELEM_WHILE_STATEMENT        = 7,
    COMPILER_PARSER_ELEM_IF_STATEMENT           = 8,
    COMPILER_PARSER_ELEM_LET_STATEMENT          = 9,
    COMPILER_PARSER_ELEM_STATEMENTS             = 10,
    COMPILER_PARSER_ELEM_SUBROUTINE_BODY        = 11,
    COMPILER_PARSER_ELEM_SUBROUTINE_DEC         = 12,
    COMPILER_PARSER_ELEM_EXPRESSION             = 13,
    COMPILER_PARSER_ELEM_TERM                   = 14,

    COMPILER_PARSER_ELEM_TOKEN,

};

char* compiler_parser_ar_elem_name[] {
    "class",                // COMPILER_PARSER_ELEM_CLASS
    "classVarDec",          // COMPILER_PARSER_ELEM_CLASS_VAR_DEC
    "parameterList",        // COMPILER_PARSER_ELEM_PARAMETER_LIST
    "varDec",               // COMPILER_PARSER_ELEM_VAR_DEC
    "expressionList",       // COMPILER_PARSER_ELEM_EXPRESSION_LIST
    "returnStatement",      // COMPILER_PARSER_ELEM_RETURN_STATEMENT
    "doStatement",          // COMPILER_PARSER_ELEM_DO_STATEMENT
    "whileStatement",       // COMPILER_PARSER_ELEM_WHILE_STATEMENT
    "ifStatement",          // COMPILER_PARSER_ELEM_IF_STATEMENT
    "letStatement",         // COMPILER_PARSER_ELEM_LET_STATEMENT
    "statements",           // COMPILER_PARSER_ELEM_STATEMENTS
    "subroutineBody",       // COMPILER_PARSER_ELEM_SUBROUTINE_BODY
    "subroutineDec",        // COMPILER_PARSER_ELEM_SUBROUTINE_DEC
    "expression",           // COMPILER_PARSER_ELEM_EXPRESSION
    "term",                 // COMPILER_PARSER_ELEM_TERM

};

struct compiler_parser_elem {
    uint8_t type;
    token_t token;

    struct compiler_parser_elem* parent;
    struct compiler_parser_elem* child;
    struct compiler_parser_elem* next;
    struct compiler_parser_elem* prev;
};

struct compiler_parser_elem* compiler_parser_add_elem(struct compiler_parser_elem* parent, uint8_t type, token_t* token) {

    struct compiler_parser_elem* cur = (struct compiler_parser_elem*)malloc(sizeof(struct compiler_parser_elem));

    if (cur == NULL)
        return NULL;

    cur->type = type;
    if (token != NULL)
        cur->token = *token;

    cur->parent = parent;
    cur->child = NULL;
    cur->prev = NULL;
    cur->next = NULL;

    if (parent != NULL) {
        if (parent->child == NULL) {
            parent->child = cur;
            cur->prev = NULL;
        } else {
            cur->prev = parent->child;
            while (cur->prev->next != NULL)
                cur->prev = cur->prev->next;
            cur->prev->next = cur;
        }
    }

    return cur;
}

void compiler_parser_print_elem(struct compiler_parser_elem* elem, int np, FILE* out_file);

void compiler_parser_print_elem(struct compiler_parser_elem* elem, int np, FILE* out_file) {
    if (elem == NULL)
        return;

    if (elem->type == COMPILER_PARSER_ELEM_TOKEN) {
        compiler_xml_print_token(&(elem->token), np, out_file);
    } else {
        compiler_xml_print_start_tag(compiler_parser_ar_elem_name[elem->type], np, out_file);
        compiler_parser_print_elem(elem->child, np+1, out_file);
        compiler_xml_print_end_tag(compiler_parser_ar_elem_name[elem->type], np, out_file);
    }

    compiler_parser_print_elem(elem->next, np, out_file);
    
    return;
}

void compiler_parser_print_elems(struct compiler_parser_elem* root, FILE* out_file) {
    compiler_parser_print_elem(root, 0, out_file);
}

int compiler_parse_is_class_var_dec(tokens_t* tokens, int* n)
{
    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL)
        return 0;

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
           (COMPILER_KEYWORD_STATIC == compiler_get_keyword_type(token->token) ||
            COMPILER_KEYWORD_FIELD == compiler_get_keyword_type(token->token))) {
            return 1;
    }

    return 0;
}

int compiler_parse_type(tokens_t* tokens, int* n, struct compiler_parser_elem* parser_elem)
{
    // type: int | char | boolean | className

    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL) {
        printf("compiler - err parser type expected type\n");
        return -1;
    }

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type) {
        if (COMPILER_KEYWORD_INT == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_CHAR == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_BOOLEAN == compiler_get_keyword_type(token->token)) {
            // nothing
         }
         else {
             printf("compiler - err parser type expected keyword int | char | boolean or className\n");
             return -1;
         }
    } else if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
        // nothing
    }
    else {
        printf("compiler - err parser type expected int | char | boolean or className\n");
        return -1;
    }

    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
    (*n)++;
    return 0;
}

int compiler_parse_class_var_dec(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // classVarDec: ('static' | 'field') type varName (',' varName)* ';'
    enum {
        state_class_var_dec_keyword           = 0,
        state_class_var_dec_type              = 1,
        state_class_var_dec_name              = 2,
        state_class_var_dec_next_name         = 3,
        state_class_var_dec_end               = 4
    }
    state = state_class_var_dec_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_CLASS_VAR_DEC, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_class_var_dec_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                        (COMPILER_KEYWORD_STATIC == compiler_get_keyword_type(token->token) ||
                         COMPILER_KEYWORD_FIELD == compiler_get_keyword_type(token->token))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_var_dec_type;
                }
                else {
                    printf("compiler - err parser class var dec expected static | field keywords\n");
                    is_break = 1;
                }
                break;

            case state_class_var_dec_type:
                if (compiler_parse_type(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_class_var_dec_name;
                break;

            case state_class_var_dec_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_var_dec_next_name;
                }
                else {
                    printf("compiler - err parser class ver dec expected name\n");
                    return -1;
                }
                break;

            case state_class_var_dec_next_name:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type) {
                    if (!strcmp(token->token, ";")) {
                        //
                        compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                        //
                        is_break = 1;
                        state = state_class_var_dec_end;
                    }
                    else if (!strcmp(token->token, ",")) {
                        //
                        compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                        //
                        state = state_class_var_dec_name;
                    }
                    else {
                        printf("compiler - err parser class ver dec expected symbols <,> or <;> after name\n");
                        return -1;
                    }
                }
                else {
                    printf("compiler - err parser class ver dec expected symbols <,> or <;> after name\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_class_var_dec_end) {
        printf("compiler: -err parser class var dec not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_param_list(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // parameterList: ((type varName) (',' (type varName)))?
    enum {
        state_param_list_empty          = 0,
        state_param_list_type           = 1,
        state_param_list_var_name       = 2,
        state_param_list_next           = 3,
        state_param_list_end            = 4
    }
    state = state_param_list_empty;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_PARAMETER_LIST, NULL); 

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_param_list_empty:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ")")) {
                    // No param
                    state = state_param_list_end;
                    is_break = 1;
                }
                else {
                    state = state_param_list_type;
                }
                (*n)--;
                break;

            case state_param_list_type:
                if (compiler_parse_type(tokens, n, parser_elem) < 0) {
                    printf("compiler: - err expected type of parameter list.\n");
                    return -1;
                }
                else {
                    state = state_param_list_var_name;
                }
                (*n)--;
                break;
                
           case state_param_list_var_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_param_list_next;
                }
                else {
                    printf("compiler - err parser parameter list expected name\n");
                    return -1;
                }
                break;
                
            case state_param_list_next:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type) {
                    if (!strcmp(token->token, ")")) {
                        state = state_param_list_end;
                        (*n)--;
                        is_break = 1;
                    }
                    else if (!strcmp(token->token, ",")) {
                        state = state_param_list_type;
                        compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    }
                    else {
                        printf("compiler - err parser type unexpected symbol %s\n", token->token);
                        return -1;
                    }
                }
                else {
                    printf("compiler - err parser parameter list expected symbol <,> or <)>\n");
                    return -1;
                }
                break;
            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_param_list_end) {
        printf("compiler: -err parser param list not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_is_var_dec(tokens_t* tokens, int* n)
{
    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL)
        return 0;

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
        COMPILER_KEYWORD_VAR == compiler_get_keyword_type(token->token)) {
            return 1;
    }

    return 0;
}

// TODO to top declaration
int compiler_parse_subroutine_call(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem);
int compiler_parse_expression(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem);

int compiler_parse_term(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // term: intConstant | stringConstant | keywordConstant | varName | '(' expression ')' |
    //       varName '[' expression ']' | (unaryOp term) | subroutineCall
    enum {
        state_term               = 0,
        state_term_expression    = 1,
        state_term_expression2   = 2,
        state_term_close_parent  = 3,
        state_term_after_name    = 4,
        state_term_term          = 5,
        state_term_close_index   = 6,
        state_term_end           = 7,
    }
    state = state_term;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_TERM, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {

            case state_term:
                // intConst | stringConstat | keywordConst
                if (COMPILER_TOKEN_NAME_INT_CONST == token->type ||
                    COMPILER_TOKEN_NAME_STRING_CONST == token->type ||
                    compiler_is_token_keyword_const(token)) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_term_end;
                }       // unaryOp
                else if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                         ((!strcmp("-", token->token)) ||
                          (!strcmp("~", token->token)))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_term_term;
                }       // '(' expression ')'
                else if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                         (!strcmp("(", token->token))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_term_expression;
                }       // var (varName | varName '[' ) or func name (funcName '(' | name '.' funcName '('
                else if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    token_t* next_token = compiler_get_token(tokens, (*n)+1);
                    if (next_token != NULL &&
                            (COMPILER_TOKEN_NAME_SYMBOL == next_token->type &&
                             ((!strcmp(next_token->token, ".")) || 
                              (!strcmp(next_token->token, "("))))) {
                        //
                        if (compiler_parse_subroutine_call(tokens, n, parser_elem) < 0)
                            return -1;
                        //
                        (*n)--;
                        is_break = 1;
                        state = state_term_end;
                    }
                    else {
                        // varName
                        compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                        //
                        state = state_term_after_name;
                    }
                }
                else {
                    printf("compiler - err parser term expected name\n");
                    return -1;
                }
                break;

            case state_term_expression:
            case state_term_expression2:
                //
                if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = (state == state_term_expression) ? state_term_close_parent : state_term_close_index;
                break;

            case state_term_term:
                //
                if (compiler_parse_term(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                is_break = 1;
                state = state_term_end;
                break;

            case state_term_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL != token->type ||
                    strcmp(")", token->token)) {
                        printf("compiler parser - term expected close ')'.\n");
                        return -1;
                }
                //
                compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                //
                is_break = 1;
                state = state_term_end;
                break;
                
            case state_term_after_name:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp("[", token->token))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token); 
                    state = state_term_expression2;
                }
                else {
                    (*n)--;
                    is_break = 1;
                    state = state_term_end;
                }
                break;

            case state_term_close_index:
                if (COMPILER_TOKEN_NAME_SYMBOL != token->type ||
                    strcmp("]", token->token)) {
                        printf("compiler parser - term expected close ']'.\n");
                        return -1;
                }
                //
                compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                //
                is_break = 1;
                state = state_term_end;
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_term_end) {
        printf("compiler: -err parser term not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_expression(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // expression: term (op term)*
    enum {
        state_expression_term          = 0,
        state_expression_op            = 1,
        state_expression_end           = 2
    }
    state = state_expression_term;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_EXPRESSION, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
           case state_expression_term:
                if (compiler_parse_term(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_expression_op;
                break;

            case state_expression_op:
                if (compiler_is_token_op(token)) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_expression_term;
                }
                else {
                    (*n)--;
                    state = state_expression_end;
                    is_break = 1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_expression_end) {
        printf("compiler: -err parser expression not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_var_dec(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // varDec: 'var' type varName (',' varName)* ';'
    enum {
        state_var_dec_keyword           = 0,
        state_var_dec_type              = 1,
        state_var_dec_name              = 2,
        state_var_dec_next              = 3,
        state_var_dec_end               = 4
    }
    state = state_var_dec_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_VAR_DEC, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_var_dec_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_VAR == compiler_get_keyword_type(token->token)) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_var_dec_type;
                }
                else {
                    printf("compiler - err parser var dec expected var keyword\n");
                    is_break = 1;
                }
                break;

           case state_var_dec_type:
                if (compiler_parse_type(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_var_dec_name;
                break;

            case state_var_dec_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_var_dec_next;
                }
                else {
                    printf("compiler - err parser var dec expected name\n");
                    is_break = 1;
                }
                break;

            case state_var_dec_next:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type) {
                    if (!strcmp(token->token, ";")) {
                        state = state_var_dec_end;
                        is_break = 1;
                    }
                    else if (!strcmp(token->token, ",")) {
                        state = state_var_dec_name;
                    }
                    else {
                        printf("compiler - err parser type unexpected symbol %s\n", token->token);
                        return -1;
                    }
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                }
                else {
                    printf("compiler - err parser type expected symbol <;> or <,>\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_var_dec_end) {
        printf("compiler: -err parser var dec not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_expression_list(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // expressionList: (expression (',' expression)*)?
    enum {
        state_expression_list_close_parent      = 0,
        state_expression_list_expression        = 1,
        state_expression_list_next              = 2,
        state_expression_list_end               = 3
    }
    state = state_expression_list_close_parent;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_EXPRESSION_LIST, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_expression_list_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ")")) {
                    //
                    
                    is_break = 1;
                    state = state_expression_list_end;
                }
                else {
                    state = state_expression_list_expression;
                }
                (*n)--;
                break;

            case state_expression_list_expression:
                if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_expression_list_next;
                break;

            case state_expression_list_next:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ",")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_expression_list_expression;
                }
                else {
                    (*n)--;
                    is_break = 1;
                    state = state_expression_list_end;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_expression_list_end) {
        printf("compiler: -err parser expressioin list not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_subroutine_call(tokens_t* tokens, int* n, struct compiler_parser_elem* parser_elem)
{
    // subroutineCall: subroutineName '(' parameterList ')' | (className | varName) '.' subroutineName '(' expressionList ')'
    enum {
        state_subroutine_call_name1             = 0,
        state_subroutine_call_dot               = 1,
        state_subroutine_call_name2             = 2,
        state_subroutine_call_start_parent      = 3,
        state_subroutine_call_expression_list   = 4,
        state_subroutine_call_close_parent      = 5,
        state_subroutine_call_end               = 7
    }
    state = state_subroutine_call_name1;

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_subroutine_call_name1:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_call_dot;
                }
                else {
                    printf("compiler - err parser subroutine call expected subroutine or calls name\n");
                    return -1;
                }
                break;

            case state_subroutine_call_dot:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ".")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_call_name2;
                }
                else {
                    (*n)--;
                    state = state_subroutine_call_start_parent;
                }
                break;

           case state_subroutine_call_name2:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_call_start_parent;
                }
                else {
                    printf("compiler - err parser subroutine call expected subroutine name\n");
                    return -1;
                }
                break;

            case state_subroutine_call_start_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "(")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_call_expression_list;
                }
                else {
                    printf("compiler - err parser subroutine call expected symbol <(>\n");
                    return -1;
                }
                break;

           case state_subroutine_call_expression_list:
                if (compiler_parse_expression_list(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_subroutine_call_close_parent;
                break;

            case state_subroutine_call_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ")")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_subroutine_call_end;
                }
                else {
                    printf("compiler - err parser subroutine call expected symbol <)>\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_subroutine_call_end) {
        printf("compiler: -err parser subroutine call not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_is_statement(tokens_t* tokens, int* n)
{
    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL)
        return 0;

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
            (COMPILER_KEYWORD_LET == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_IF == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_WHILE == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_DO == compiler_get_keyword_type(token->token) ||
             COMPILER_KEYWORD_RETURN == compiler_get_keyword_type(token->token))) {
            return 1;
    }

    return 0;
}


int compiler_parse_return_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // returnStatement: 'return' expression? ';'
    enum {
        state_return_statement_keyword          = 0,
        state_return_statement_expression       = 1,
        state_return_statement_terminator       = 2,
        state_return_statement_end              = 3
    }
    state = state_return_statement_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_RETURN_STATEMENT, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_return_statement_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_RETURN == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    state = state_return_statement_expression;
                }
                else {
                    return -1;
                }
                break;

            case state_return_statement_expression:
               if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ";"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_return_statement_end;
                }
                else {
                    if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                        return -1;
                    }
                    (*n)--;
                    state = state_return_statement_terminator;
                }
                break;
                
            case state_return_statement_terminator:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ";"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_return_statement_end;
                }

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_return_statement_end) {
        printf("compiler: -err parser return statement not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_do_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // doStatement: 'do' subroutineCall ';'
    enum {
        state_do_statement_keyword          = 0,
        state_do_statement_subroutine_call  = 1,
        state_do_statement_terminator       = 2,
        state_do_statement_end              = 3
    }
    state = state_do_statement_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_DO_STATEMENT, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_do_statement_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_DO == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    state = state_do_statement_subroutine_call;
                }
                else {
                    return -1;
                }
                break;

            case state_do_statement_subroutine_call:
                if (compiler_parse_subroutine_call(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_do_statement_terminator;
                break;

            case state_do_statement_terminator:
               if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ";"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_do_statement_end;
                }
                else {
                    printf("compiler - err parser let statement expected <;> symbol\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_do_statement_end) {
        printf("compiler: -err parser do statement not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

// TODO - in top declaration
int compiler_parse_statements(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem);

int compiler_parse_while_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // whileStatement: 'while' '(' expression ')' '{' statements '}'
    enum {
        state_while_statement_keyword          = 0,
        state_while_statement_open_parent      = 1,
        state_while_statement_expression       = 2,
        state_while_statement_close_parent     = 3,
        state_while_statement_open_body        = 4,
        state_while_statement_statements       = 5,
        state_while_statement_close_body       = 6,
        state_while_statement_end              = 11
    }
    state = state_while_statement_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_WHILE_STATEMENT, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_while_statement_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_WHILE == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    state = state_while_statement_open_parent;
                }
                else {
                    return -1;
                }
                break;

            case state_while_statement_open_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "("))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_while_statement_expression;
                }
                else {
                    printf("compiler - err parser while statement expected <(> symbol\n");
                    return -1;
                }
                break;

            case state_while_statement_expression:
                if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_while_statement_close_parent;
                break;

            case state_while_statement_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ")"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_while_statement_open_body;
                }
                else {
                    printf("compiler - err parser while statement expected <)> symbol\n");
                    return -1;
                }
                break;

            case state_while_statement_open_body:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "{"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_while_statement_statements;
                }
                else {
                    printf("compiler - err parser while statement expected <{> symbol\n");
                    return -1;
                }
                break;

           case state_while_statement_statements:
                if (compiler_parse_statements(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                state = state_while_statement_close_body;
                (*n)--;
                break;

            case state_while_statement_close_body:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "}"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_while_statement_end;
                }
                else {
                    printf("compiler - err parser while statement expected <}> symbol\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_while_statement_end) {
        printf("compiler: -err parser while statement not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_if_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // ifStatement: 'if' '(' expression ')' '{' statements '}' ('else' '{' statements '}')?
    enum {
        state_if_statement_keyword          = 0,
        state_if_statement_open_parent      = 1,
        state_if_statement_expression       = 2,
        state_if_statement_close_parent     = 3,
        state_if_statement_open_body        = 4,
        state_if_statement_statements       = 5,
        state_if_statement_close_body       = 6,
        state_if_statement_else             = 7,
        state_if_statement_end              = 8
    }
    state = state_if_statement_keyword;
    uint8_t is_else = 0;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_IF_STATEMENT, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_if_statement_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_IF == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    state = state_if_statement_open_parent;
                }
                else {
                    return -1;
                }
                break;

            case state_if_statement_open_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "("))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_if_statement_expression;
                }
                else {
                    printf("compiler - err parser if statement expected <(> symbol\n");
                    return -1;
                }
                break;

            case state_if_statement_expression:
                if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = state_if_statement_close_parent;
                break;

            case state_if_statement_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ")"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_if_statement_open_body;
                }
                else {
                    printf("compiler - err parser if statement expected <)> symbol\n");
                    return -1;
                }
                break;

            case state_if_statement_open_body:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "{"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_if_statement_statements;
                }
                else {
                    printf("compiler - err parser if statement expected <{> symbol\n");
                    return -1;
                }
                break;

           case state_if_statement_statements:
                if (compiler_parse_statements(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                state = state_if_statement_close_body;
                (*n)--;
                break;

            case state_if_statement_close_body:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "}"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    if (is_else) {
                        is_break = 1;
                        state = state_if_statement_end;
                    }
                    else {
                        state = state_if_statement_else;
                    }
                }
                else {
                    printf("compiler - err parser if statement expected <}> symbol\n");
                    return -1;
                }
                break;

            case state_if_statement_else:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_ELSE == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    is_else = 1;
                    state = state_if_statement_open_body;
                }
                else {
                    (*n)--;
                    is_break = 1;
                    state = state_if_statement_end;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_if_statement_end) {
        printf("compiler: -err parser if statement not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_let_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // letStatement: 'let' varName ( '[' expression ']' )? '=' expression ';'
    enum {
        state_let_statement_keyword          = 0,
        state_let_statement_var_name         = 1,
        state_let_statement_after_var_name   = 2,
        state_let_statement_index_expression = 3,
        state_let_statement_close_index      = 4,
        state_let_statement_assignment      = 5,
        state_let_statement_expression       = 6,
        state_let_statement_terminator       = 7,
        state_let_statement_end              = 8
    }
    state = state_let_statement_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_LET_STATEMENT, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_let_statement_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_LET == compiler_get_keyword_type(token->token)) {
                    
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    state = state_let_statement_var_name;
                }
                else {
                    return -1;
                }
                break;

            case state_let_statement_var_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_let_statement_after_var_name;
                }
                else {
                    printf("compiler - err parser let statement expected var name\n");
                    return -1;
                }
                break;

            case state_let_statement_after_var_name:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "["))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_let_statement_index_expression;
                }
                else {
                    (*n)--;
                    state = state_let_statement_assignment;
                }
                break;

            case state_let_statement_assignment:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "="))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_let_statement_expression;
                }
                else {
                    printf("compiler - err parser let statement expected <=> symbol\n");
                    return -1;
                }
                break;

            case state_let_statement_expression:
            case state_let_statement_index_expression:
                if (compiler_parse_expression(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                state = (state == state_let_statement_expression) ? 
                            state_let_statement_terminator : state_let_statement_close_index;
                break;

            case state_let_statement_close_index:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, "]"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_let_statement_assignment;
                }
                else {
                    printf("compiler - err parser let statement expected <]> symbol\n");
                    return -1;
                }
                break;

            case state_let_statement_terminator:
               if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    (!strcmp(token->token, ";"))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    is_break = 1;
                    state = state_let_statement_end;
                }
                else {
                    printf("compiler - err parser let statement expected <;> symbol\n");
                    return -1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_let_statement_end) {
        printf("compiler: -err parser let statement not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_statement(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // statement: letStatement | ifStatement | whileStatement | doStatement | returnStatement

    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL) {
        printf("compiler - err parser statement expected keyword\n");
        return -1;
    }

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type) {
        if (COMPILER_KEYWORD_LET == compiler_get_keyword_type(token->token)) {
            compiler_parse_let_statement(tokens, n, parent_elem);
        }
        else if (COMPILER_KEYWORD_IF == compiler_get_keyword_type(token->token)) {
            compiler_parse_if_statement(tokens, n, parent_elem);
        }
        else if (COMPILER_KEYWORD_WHILE == compiler_get_keyword_type(token->token)) {
            compiler_parse_while_statement(tokens, n, parent_elem);
        }
        else if (COMPILER_KEYWORD_DO == compiler_get_keyword_type(token->token)) {
            compiler_parse_do_statement(tokens, n, parent_elem);
        }
        else if (COMPILER_KEYWORD_RETURN == compiler_get_keyword_type(token->token)) {
            compiler_parse_return_statement(tokens, n, parent_elem);
        }
        else {
            printf("compiler - err parser statement expected let | if | while | do | return keyword\n");
            return -1;
        }
    }
    else {
        printf("compiler - err parser statement expected keyword\n");
        return -1;
    }

    return 0;
}

int compiler_parse_statements(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // statements: statement*
    enum {
        state_statements_statement         = 0,
        state_statements_end               = 1
    }
    state = state_statements_statement;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_STATEMENTS, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_statements_statement:
                if (compiler_parse_is_statement(tokens, n)) {
                    if (compiler_parse_statement(tokens, n, parser_elem) < 0) {
                        return -1;
                    }
                    (*n)--;
                }
                else {
                    (*n)--;
                    is_break = 1;
                    state = state_statements_end;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_statements_end) {
        printf("compiler: -err parser statements not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}


int compiler_parse_subroutine_body(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // subroutineBody: '{' varDec* statements '}'
    enum {
        state_subroutine_body_open_symbol  = 0,
        state_subroutine_body_var_dec      = 1,
        state_subroutine_body_statements   = 2,
        state_subroutine_body_close_symbol = 3,
        state_subroutine_body_end          = 4
    }
    state = state_subroutine_body_open_symbol;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_SUBROUTINE_BODY, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_subroutine_body_open_symbol:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "{")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_body_var_dec;
                }
                else {
                    printf("compiler - err parser subroutine body expected symbol <{>\n");
                    return -1;
                }
                break;

           case state_subroutine_body_var_dec:
                if (compiler_parse_is_var_dec(tokens, n)) {
                    if (compiler_parse_var_dec(tokens, n, parser_elem) < 0) {
                        return -1;
                    }
                }
                else {
                    state = state_subroutine_body_statements;
                }
                (*n)--;
                break;

           case state_subroutine_body_statements:
                if (compiler_parse_statements(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                state = state_subroutine_body_close_symbol;
                (*n)--;
                break;

            case state_subroutine_body_close_symbol:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "}")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_body_end;
                    is_break = 1;
                }
                else {
                    printf("compiler - err parser subroutine body expected symbol <}>\n");
                    return -1;
                }
                break;
                
                // TODO

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_subroutine_body_end) {
        printf("compiler: -err parser subroutine not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_is_subroutine_dec(tokens_t* tokens, int* n)
{
    token_t* token = compiler_get_token(tokens, *n);

    if (token == NULL)
        return 0;

    if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
           (COMPILER_KEYWORD_CONSTRUCTOR == compiler_get_keyword_type(token->token) ||
            COMPILER_KEYWORD_METHOD == compiler_get_keyword_type(token->token) ||
            COMPILER_KEYWORD_FUNCTION == compiler_get_keyword_type(token->token))) {
            return 1;
    }

    return 0;
}

int compiler_parse_subroutine_dec(tokens_t* tokens, int* n, struct compiler_parser_elem* parent_elem)
{
    // subroutineDec: ('constructor' | 'function' | 'method') ('void' | type) 
    //                  subroutineName '(' parameterList ')' subroutineBody
    enum {
        state_subroutine_dec_keyword           = 0,
        state_subroutine_dec_type              = 1,
        state_subroutine_dec_name              = 2,
        state_subroutine_dec_start_parent      = 3,
        state_subroutine_dec_param_list        = 4,
        state_subroutine_dec_close_parent      = 5,
        state_subroutine_dec_body              = 6,
        state_subroutine_dec_end               = 7
    }
    state = state_subroutine_dec_keyword;

    struct compiler_parser_elem* parser_elem = compiler_parser_add_elem(parent_elem, COMPILER_PARSER_ELEM_SUBROUTINE_DEC, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_subroutine_dec_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                        (COMPILER_KEYWORD_CONSTRUCTOR == compiler_get_keyword_type(token->token) ||
                         COMPILER_KEYWORD_METHOD == compiler_get_keyword_type(token->token) ||
                         COMPILER_KEYWORD_FUNCTION == compiler_get_keyword_type(token->token))) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_dec_type;
                }
                else {
                    printf("compiler - err parser subroutine dec expected constructor | function | method keywords\n");
                    is_break = 1;
                }
                break;

            case state_subroutine_dec_type:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                        COMPILER_KEYWORD_VOID == compiler_get_keyword_type(token->token)) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                }
                else {
                    if (compiler_parse_type(tokens, n, parser_elem) < 0) {
                        return -1;
                    }
                    (*n)--;
                }
                state = state_subroutine_dec_name;
                break;

            case state_subroutine_dec_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_dec_start_parent;
                }
                else {
                    printf("compiler - err parser subroutine dec expected name\n");
                    is_break = 1;
                }
                break;

            case state_subroutine_dec_start_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "(")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_dec_param_list;
                }
                else {
                    printf("compiler - err parser subroutine dec expected symbol <(>\n");
                    is_break = 1;
                }
                break;

           case state_subroutine_dec_param_list:
                if (compiler_parse_param_list(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                state = state_subroutine_dec_close_parent;
                (*n)--;
                break;

            case state_subroutine_dec_close_parent:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, ")")) {
                    //
                    compiler_parser_add_elem(parser_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_subroutine_dec_body;
                }
                else {
                    printf("compiler - err parser subroutine dec expected symbol <)>\n");
                    is_break = 1;
                }
                break;

           case state_subroutine_dec_body:
                if (compiler_parse_subroutine_body(tokens, n, parser_elem) < 0) {
                    return -1;
                }
                (*n)--;
                is_break = 1;
                state = state_subroutine_dec_end;
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_subroutine_dec_end) {
        printf("compiler: -err parser subroutineDec not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

int compiler_parse_class(tokens_t* tokens, int* n, struct compiler_parser_elem** pclass_elem)
{
    // class: 'class' className '{' classVarDec* subroutineDec* '}'
    enum {
        state_class_keyword           = 0,
        state_class_name              = 1,
        state_class_body_open_symbol  = 2,
        state_class_var_dec           = 3,
        state_class_subroutine_dec    = 4,
        state_class_body_close_symbol = 5,
        state_class_end               = 6
    }
    state = state_class_keyword;

    *pclass_elem = compiler_parser_add_elem(NULL, COMPILER_PARSER_ELEM_CLASS, NULL);

    token_t* token;
    while ((token = compiler_get_token(tokens, *n)) != NULL) {
        uint8_t is_break = 0;

        switch (state) {
            case state_class_keyword:
                if (COMPILER_TOKEN_NAME_KEYWORD == token->type &&
                    COMPILER_KEYWORD_CLASS == compiler_get_keyword_type(token->token)) {
                    //
                    compiler_parser_add_elem(*pclass_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_name;
                }
                else {
                    printf("compiler - err parser class expected class keywords\n");
                    is_break = 1;
                }
                break;

            case state_class_name:
                if (COMPILER_TOKEN_NAME_IDENTIFIER == token->type) {
                    //
                    compiler_parser_add_elem(*pclass_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_body_open_symbol;
                }
                else {
                    printf("compiler - err parser class expected class name\n");
                    is_break = 1;
                }
                break;

            case state_class_body_open_symbol:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "{")) {
                    //
                    compiler_parser_add_elem(*pclass_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_var_dec;
                }
                else {
                    printf("compiler - err parser class expected symbol <{>\n");
                    is_break = 1;
                }
                break;

            case state_class_var_dec:
                if (compiler_parse_is_class_var_dec(tokens, n)) {
                    if (compiler_parse_class_var_dec(tokens, n, *pclass_elem) < 0) {
                        return -1;
                    }
                }
                else {
                    state = state_class_subroutine_dec;
                }
                (*n)--;
                break;

            case state_class_subroutine_dec:
                if (compiler_parse_is_subroutine_dec(tokens, n)) {
                    if (compiler_parse_subroutine_dec(tokens, n, *pclass_elem) < 0) {
                        return -1;
                    }
                }
                else {
                    state = state_class_body_close_symbol;
                }
                (*n)--;
                break;

            case state_class_body_close_symbol:
                if (COMPILER_TOKEN_NAME_SYMBOL == token->type &&
                    !strcmp(token->token, "}")) {
                    //
                    compiler_parser_add_elem(*pclass_elem, COMPILER_PARSER_ELEM_TOKEN, token);
                    //
                    state = state_class_end;
                }
                else {
                    printf("compiler - err parser class expected symbol <}>\n");
                    is_break = 1;
                }
                break;

            default:
                break;
        }

        (*n)++;
        if (is_break)
            break;
    }

    if (state != state_class_end) {
        printf("compiler: -err parser class not finish, state = %d\n", state);
        return -1;
    }

    return 0;
}

token_t* compiler_code_gen_get_class_name_token(struct compiler_parser_elem* class_elem)
{
    if (class_elem == NULL || class_elem->type != COMPILER_PARSER_ELEM_CLASS) {
        printf("compiler code gen - error absent class\n");
        return NULL;
    }

    struct compiler_parser_elem* parse_elem = class_elem->child;

    if (parse_elem == NULL) {
        printf("compiler code gen - error empty class parse element\n");
        return NULL;
    }

    if (parse_elem->type != COMPILER_PARSER_ELEM_TOKEN ||
        parse_elem->token.type != COMPILER_TOKEN_NAME_KEYWORD ||
        COMPILER_KEYWORD_CLASS != compiler_get_keyword_type(parse_elem->token.token)) {
        printf("compiler code gen - error class parse element without 'class' keyword\n");
        return NULL;
    }

    parse_elem = parse_elem->next;
    if (parse_elem == NULL ||
        parse_elem->token.type != COMPILER_TOKEN_NAME_IDENTIFIER) {
        printf("compiler code gen - error class parse element without className\n");
        return NULL;
    }

    return &(parse_elem->token);
}

enum compiler_code_gen_routine_kind {
    COMPILER_CODE_GEN_ROUTINE_KIND_CONSTRUCTOR  = 0,
    COMPILER_CODE_GEN_ROUTINE_KIND_METHOD       = 1,
    COMPILER_CODE_GEN_ROUTINE_KIND_FUNCTION     = 2,
    COMPILER_CODE_GEN_ROUTINE_KIND_ERROR        = 3
};

char* compiler_code_gen_routine_ar_kind_name[] = {
    "constructor",        // = COMPILER_CODE_GEN_ROUTINE_KIND_CONSTRUCTOR
    "method",             // = COMPILER_CODE_GEN_ROUTINE_KIND_METHOD
    "function",           // = COMPILER_CODE_GEN_ROUTINE_KIND_FUNCTION
    "_error_",            // = COMPILER_CODE_GEN_ROUTINE_KIND_ERROR
};

uint8_t compiler_code_gen_get_routine_kind(struct compiler_parser_elem* fun_elem)
{
    uint8_t routine_kind = COMPILER_CODE_GEN_ROUTINE_KIND_ERROR;
    
    if (fun_elem == NULL ||
        fun_elem->type != COMPILER_PARSER_ELEM_SUBROUTINE_DEC) {
        return routine_kind;
    }

    if (fun_elem->child != NULL /* kind */ &&
        COMPILER_TOKEN_NAME_KEYWORD == fun_elem->child->token.type) {
        //
        uint8_t type = compiler_get_keyword_type(fun_elem->child->token.token);
        switch (type) {
            case COMPILER_KEYWORD_CONSTRUCTOR:
                routine_kind = COMPILER_CODE_GEN_ROUTINE_KIND_CONSTRUCTOR;
                break;
            case COMPILER_KEYWORD_METHOD:
                routine_kind = COMPILER_CODE_GEN_ROUTINE_KIND_METHOD;
                break;
            case COMPILER_KEYWORD_FUNCTION:
                routine_kind = COMPILER_CODE_GEN_ROUTINE_KIND_FUNCTION;
                break;
        }
    }

    return routine_kind;
}

token_t* compiler_code_gen_get_next_fun(struct compiler_parser_elem** p_fun_elem)
{
    if ((*p_fun_elem) == NULL) {
        return NULL;
    }

    if ((*p_fun_elem)->type == COMPILER_PARSER_ELEM_CLASS) {
        // First in class
        (*p_fun_elem) = (*p_fun_elem)->child;
    }
    else {
        // Next function.
        *p_fun_elem = (*p_fun_elem)->next;
    }

    // Find a function
    while (*p_fun_elem != NULL) {
        if ((*p_fun_elem)->type == COMPILER_PARSER_ELEM_SUBROUTINE_DEC)
            break;
        *p_fun_elem = (*p_fun_elem)->next;
    }
    if (*p_fun_elem == NULL)
        return NULL;

    if ((*p_fun_elem)->child != NULL &&             /* kind */ 
        (*p_fun_elem)->child->next != NULL &&       /* type */
        (*p_fun_elem)->child->next->next != NULL) { /* name */
           return &((*p_fun_elem)->child->next->next->token);
    }

    printf("compiler code gen - err function without name.\n");
    return NULL;
}

enum compiler_code_gen_var_kind {
    COMPILER_CODE_GEN_VAR_KIND_FIELD        = 0,
    COMPILER_CODE_GEN_VAR_KIND_STATIC       = 1,
    COMPILER_CODE_GEN_VAR_KIND_ARG          = 2,
    COMPILER_CODE_GEN_VAR_KIND_VAR          = 3,
};

char* compiler_code_gen_var_ar_kind_name[] = {
    "field",        // = COMPILER_CODE_GEN_VAR_KIND_FIELD
    "static",       // = COMPILER_CODE_GEN_VAR_KIND_STATIC
    "arg",          // = COMPILER_CODE_GEN_VAR_KIND_ARG
    "var",          // = COMPILER_CODE_GEN_VAR_KIND_VAR
};

struct compiler_code_gen_var {
    char*       name;
    char*       type;
    uint8_t     kind;
    int         num;
};

#define COMPILER_CODE_GEN_VAR_TABLE_MAX_NUM 1024

enum compiler_code_gen_var_section {
    COMPILER_CODE_GEN_VAR_SECTION_THIS        = 0,
    COMPILER_CODE_GEN_VAR_SECTION_STATIC      = 1,
    COMPILER_CODE_GEN_VAR_SECTION_ARGUMENT    = 2,
    COMPILER_CODE_GEN_VAR_SECTION_LOCAL       = 3,

    COMPILER_CODE_GEN_VAR_SECTION_UNKNOWN     = 4,
};

char* compiler_code_gen_var_ar_section_name[] = {
    "this",        // = COMPILER_CODE_GEN_VAR_SECTION_THIS
    "static",      // = COMPILER_CODE_GEN_VAR_SECTION_STATIC
    "argument",    // = COMPILER_CODE_GEN_VAR_SECTION_ARGUMENT
    "local",       // = COMPILER_CODE_GEN_VAR_SECTION_LOCAL

    "unknown",     // = COMPILER_CODE_GEN_VAR_SECTION_UNKNOWN
};

uint8_t compiler_code_gen_get_var_section(struct compiler_code_gen_var* var) {
    uint8_t section = COMPILER_CODE_GEN_VAR_SECTION_UNKNOWN;
    switch (var->kind) {
        case COMPILER_CODE_GEN_VAR_KIND_FIELD:
            section = COMPILER_CODE_GEN_VAR_SECTION_THIS;
            break;
        case COMPILER_CODE_GEN_VAR_KIND_STATIC:
            section = COMPILER_CODE_GEN_VAR_SECTION_STATIC;
            break;
        case COMPILER_CODE_GEN_VAR_KIND_ARG:
            section = COMPILER_CODE_GEN_VAR_SECTION_ARGUMENT;
            break;
        case COMPILER_CODE_GEN_VAR_KIND_VAR:
            section = COMPILER_CODE_GEN_VAR_SECTION_LOCAL;
            break;
    
    }
    return section;
}

struct compiler_code_var_table {
    struct compiler_code_gen_var symbols[COMPILER_CODE_GEN_VAR_TABLE_MAX_NUM];
    int num;
} compiler_code_var_class_table, compiler_code_var_subroutine_table;

void compiler_code_gen_init_var_table(struct compiler_code_var_table* table) {
    table->num = 0;
}

int compiler_code_gen_add_to_var_table(struct compiler_code_var_table* table, struct compiler_code_gen_var* var) {
    if (table->num >= COMPILER_CODE_GEN_VAR_TABLE_MAX_NUM)
        return -1;

    table->symbols[(table->num)++] = *var;
    return 0;
}

struct compiler_code_gen_var* compiler_code_gen_get_from_var_table(struct compiler_code_var_table* table, int n) {
    if (n >= table->num)
        return NULL;

    return &(table->symbols[n]);
}

struct compiler_code_gen_var* compiler_code_gen_find_from_var_table(struct compiler_code_var_table* table,  char* var_name) {
    
    int i = 0;
    while (i < table->num) {
        if (!strcmp(var_name, table->symbols[i].name))
            break;
        ++i;
    }

    if (i == table->num)
        return NULL;

    return &(table->symbols[i]);
}

struct compiler_code_gen_var* compiler_code_gen_find_fun_var(struct compiler_code_var_table* class_table,
                                                             struct compiler_code_var_table* fun_table,
                                                             char* var_name) {

    struct compiler_code_gen_var* var = compiler_code_gen_find_from_var_table(fun_table, var_name);

    if (var != NULL)
        return var;

    return compiler_code_gen_find_from_var_table(class_table, var_name);
}

int compiler_code_gen_vm_var(char* class_name,
                             struct compiler_code_gen_var* var,
                             FILE* out_file)
{
    char print_line[2048+1];
    uint8_t section = compiler_code_gen_get_var_section(var);

    //
    sprintf(print_line, " %s",
                compiler_code_gen_var_ar_section_name[section]);
    fwrite(print_line, 1, strlen(print_line), out_file);

    if (section == COMPILER_CODE_GEN_VAR_SECTION_STATIC) {
        sprintf(print_line, " %s.%s\n",
                    class_name, var->name);
    } else {
        sprintf(print_line, " %d\n", var->num);
    }
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

int compiler_code_gen_vm_pop_var(char* class_name,
                                 struct compiler_code_gen_var* var,
                                 FILE* out_file)
{
    char print_line[2048+1];
    //
    sprintf(print_line, "pop");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return compiler_code_gen_vm_var(class_name, var, out_file);
}

int compiler_code_gen_vm_push_var(char* class_name,
                                  struct compiler_code_gen_var* var,
                                  FILE* out_file)
{
    char print_line[2048+1];
    //
    sprintf(print_line, "push");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return compiler_code_gen_vm_var(class_name, var, out_file);
}

int compiler_gen_code_parse_class_var( struct compiler_parser_elem** p_parse_elem,
                                       struct compiler_code_var_table* table, 
                                       int* n_field, int* n_static)
{
    if (*p_parse_elem == NULL)
        return 0;

    if ((*p_parse_elem)->type == COMPILER_PARSER_ELEM_CLASS)
        *p_parse_elem = (*p_parse_elem)->child;

    while (*p_parse_elem != NULL) {
        if ((*p_parse_elem)->type == COMPILER_PARSER_ELEM_CLASS_VAR_DEC)
            break;
        *p_parse_elem = (*p_parse_elem)->next;
    }
    if (*p_parse_elem == NULL)
        return 0;

    int n = 0;
    struct compiler_code_gen_var var;
    struct compiler_parser_elem* var_elem = (*p_parse_elem)->child;
    if (var_elem == NULL) {
        printf("compiler code gen - err empty class var_elem\n");
        return -1;
    }
    // Kind
    if (var_elem->token.type != COMPILER_TOKEN_NAME_KEYWORD ||
        (COMPILER_KEYWORD_FIELD  != compiler_get_keyword_type(var_elem->token.token) &&
         COMPILER_KEYWORD_STATIC != compiler_get_keyword_type(var_elem->token.token))) {
        printf("compiler code gen - err kind = %s of class var_elem\n", var_elem->token.token);
        return -1;
    }

    var.kind = (COMPILER_KEYWORD_FIELD == compiler_get_keyword_type(var_elem->token.token)) ?
                        COMPILER_CODE_GEN_VAR_KIND_FIELD : COMPILER_CODE_GEN_VAR_KIND_STATIC;

    // Type
    if ((var_elem = var_elem->next) == NULL) {
        printf("compiler code gen - err type of class var_elem\n");
        return -1;
    }
    var.type = var_elem->token.token;

    while (1) {
        if (n > 0) {
            // check ';' or ','
            if ((var_elem = var_elem->next) == NULL ||
                var_elem->token.type != COMPILER_TOKEN_NAME_SYMBOL ||
                ((!strcmp(",", var_elem->token.token)) &&
                 (!strcmp(";", var_elem->token.token)))) {
                printf("compiler code gen - err expected ',' or ';' after name.\n");
                return -1;
            }
            if (!strcmp(";", var_elem->token.token))
               break;
        }

        // Name
        if ((var_elem = var_elem->next) == NULL ||
            var_elem->token.type != COMPILER_TOKEN_NAME_IDENTIFIER) {
            printf("compiler code gen - err name of class var_elem\n");
            return -1;
        }
        var.name = var_elem->token.token;

        // Num
        var.num = (var.kind == COMPILER_CODE_GEN_VAR_KIND_FIELD) ? (*n_field)++ : (*n_static)++;

        // Add a var
        printf("=== %d Add var.name = %s\n", n, var.name);
        if (compiler_code_gen_add_to_var_table(table, &var) < 0) {
            printf("compiler: code gen - error to add var to class symbol table\n");
            return -1;
        }
        ++n;
    }

    *p_parse_elem = (*p_parse_elem)->next;
    return n;
}

int compiler_gen_code_parse_fun_arg_var(struct compiler_parser_elem* parse_elem,
                                        struct compiler_code_var_table* table)
{
    if (parse_elem == NULL ||
        parse_elem->type != COMPILER_PARSER_ELEM_SUBROUTINE_DEC) {
        printf("compiler code gen - err find arg vars not in subroutine.\n");
        return -1;
    }

    parse_elem = parse_elem->child;

    while (parse_elem != NULL) {
        if (parse_elem->type == COMPILER_PARSER_ELEM_PARAMETER_LIST)
            break;
        parse_elem = parse_elem->next;
    }
    if (parse_elem == NULL) {
        printf("compiler code gen - arg vars absent parameter list.\n");
        return -1;
    }

    int n = 0;
    struct compiler_code_gen_var var;
    parse_elem = parse_elem->child;
    if (parse_elem == NULL)
        return 0;

    // Kind
    var.kind = COMPILER_CODE_GEN_VAR_KIND_ARG;

    while (1) {
        if (n > 0) {
            // check ','
            if ((parse_elem = parse_elem->next) == NULL)
                break;
            if (parse_elem->token.type != COMPILER_TOKEN_NAME_SYMBOL ||
                strcmp(",", parse_elem->token.token)) {
                printf("compiler code gen - err expected ',' after name.\n");
                return -1;
            }
        }

        // Type
        if (n > 0)
            parse_elem = parse_elem->next;
        var.type = parse_elem->token.token;

        // Name
        if ((parse_elem = parse_elem->next) == NULL ||
            parse_elem->token.type != COMPILER_TOKEN_NAME_IDENTIFIER) {
            printf("compiler code gen - err name of parameter var.\n");
            return -1;
        }
        var.name = parse_elem->token.token;

        // Num
        var.num = n;

        // Add a var
        printf("=== %d Add arg var.name = %s\n", n, var.name);
        if (compiler_code_gen_add_to_var_table(table, &var) < 0) {
            printf("compiler: code gen - error to add arg var to fun symbol table\n");
            return -1;
        }
        ++n;
    }

    return n;
}

int compiler_gen_code_parse_fun_var(struct compiler_parser_elem** p_parse_elem,
                                    struct compiler_code_var_table* table, 
                                    int* n_var)
{
    if (*p_parse_elem == NULL)
        return 0;

    if ((*p_parse_elem)->type == COMPILER_PARSER_ELEM_SUBROUTINE_DEC) {
        //
        *p_parse_elem = (*p_parse_elem)->child;

        while (*p_parse_elem != NULL) {
            if ((*p_parse_elem)->type == COMPILER_PARSER_ELEM_SUBROUTINE_BODY)
               break;
            *p_parse_elem = (*p_parse_elem)->next;
        }
        if (*p_parse_elem == NULL)
            return 0;

        *p_parse_elem = (*p_parse_elem)->child;
    }

    while (*p_parse_elem != NULL) {
        if ((*p_parse_elem)->type == COMPILER_PARSER_ELEM_VAR_DEC)
            break;
        *p_parse_elem = (*p_parse_elem)->next;
    }
    if (*p_parse_elem == NULL)
        return 0;

    int n = 0;
    struct compiler_code_gen_var var;
    struct compiler_parser_elem* var_elem = (*p_parse_elem)->child;
    if (var_elem == NULL) {
        printf("compiler code gen - err var dec is empty.\n");
        return -1;
    }

    // Kind
    if (var_elem->token.type != COMPILER_TOKEN_NAME_KEYWORD ||
        COMPILER_KEYWORD_VAR  != compiler_get_keyword_type(var_elem->token.token)) {
        printf("compiler code gen - err find <%s> instead 'var' in var dec.\n", var_elem->token.token);
        return -1;
    }

    var.kind = COMPILER_CODE_GEN_VAR_KIND_VAR;

    // Type
    if ((var_elem = var_elem->next) == NULL) {
        printf("compiler code gen - err type of var dec.\n");
        return -1;
    }
    var.type = var_elem->token.token;

    while (1) {
        if (n > 0) {
            // check ';' or ','
            if ((var_elem = var_elem->next) == NULL ||
                var_elem->token.type != COMPILER_TOKEN_NAME_SYMBOL ||
                ((!strcmp(",", var_elem->token.token)) &&
                 (!strcmp(";", var_elem->token.token)))) {
                printf("compiler code gen - err expected ',' or ';' after name.\n");
                return -1;
            }
            if (!strcmp(";", var_elem->token.token))
               break;
        }

        // Name
        if ((var_elem = var_elem->next) == NULL ||
            var_elem->token.type != COMPILER_TOKEN_NAME_IDENTIFIER) {
            printf("compiler code gen - err name of class var_elem\n");
            return -1;
        }
        var.name = var_elem->token.token;

        // Num
        var.num = (*n_var)++;

        // Add a var
        printf("=== %d Add var.name = %s\n", n, var.name);
        if (compiler_code_gen_add_to_var_table(table, &var) < 0) {
            printf("compiler: code gen - error to add var to function symbol table\n");
            return -1;
        }
        ++n;
    }

    *p_parse_elem = (*p_parse_elem)->next;
    return n;
}

// TODO to top
int compiler_code_gen_vm_expression(char* class_name,
                                    char* fun_name,
                                    uint8_t routine_kind,
                                    struct compiler_parser_elem* expression_elem,
                                    struct compiler_code_var_table* var_class_table,
                                    struct compiler_code_var_table* var_subroutine_table,
                                    int n_class_field,
                                    int n_fun_var,
                                    FILE* out_file);

int compiler_code_gen_vm_term(char* class_name,
                              char* fun_name,
                              uint8_t routine_kind,
                              struct compiler_parser_elem* term_elem,
                              struct compiler_code_var_table* var_class_table,
                              struct compiler_code_var_table* var_subroutine_table,
                              int n_class_field,
                              int n_fun_var,
                              FILE* out_file)
{
    // term: intConstant | stringConstant | keywordConstant | varName | '(' expression ')' |
    //       varName '[' expression ']' | (unaryOp term) | subroutineCall

    char print_line[2048+1];

    //
    sprintf(print_line, "// term\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

int compiler_code_gen_vm_expression(char* class_name,
                                    char* fun_name,
                                    uint8_t routine_kind,
                                    struct compiler_parser_elem* expression_elem,
                                    struct compiler_code_var_table* var_class_table,
                                    struct compiler_code_var_table* var_subroutine_table,
                                    int n_class_field,
                                    int n_fun_var,
                                    FILE* out_file)
{
    // expression: term (op term)*

    char print_line[2048+1];

    //
    sprintf(print_line, "// expression\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

// TODO to top
int compiler_code_gen_vm_statements(char* class_name,
                                    char* fun_name,
                                    uint8_t routine_kind,
                                    struct compiler_parser_elem* statements_elem,
                                    struct compiler_code_var_table* var_class_table,
                                    struct compiler_code_var_table* var_subroutine_table,
                                    int n_class_field,
                                    int n_fun_var,
                                    FILE* out_file);

int compiler_code_gen_vm_let_statement(char* class_name,
                                       char* fun_name,
                                       uint8_t routine_kind,
                                       struct compiler_parser_elem* statement_elem,
                                       struct compiler_code_var_table* var_class_table,
                                       struct compiler_code_var_table* var_subroutine_table,
                                       int n_class_field,
                                       int n_fun_var,
                                       FILE* out_file)
{
    // letStatement: 'let' varName ( '[' expression ']' )? '=' expression ';'

    char print_line[2048+1];

    //
    sprintf(print_line, "// let statement\n");
    fwrite(print_line, 1, strlen(print_line), out_file);


    if (statement_elem == NULL) {
        printf("Compiler code gen vm let statement - err NULL  parse elem.\n");
        return -1;
    }
    if (statement_elem->type != COMPILER_PARSER_ELEM_LET_STATEMENT) {
        printf("Compiler code gen vm let statement - err not let statement in parse elem.\n");
        return -1;
    }
    statement_elem = statement_elem->child;

    if (statement_elem == NULL ||
        statement_elem->next == NULL ||
        statement_elem->next->next == NULL ||
        statement_elem->next->next->next == NULL) {
        printf("Compiler code gen vm let statement - err < 4 elems.\n");
        return -1;
    }

    if (statement_elem->next->next->type == COMPILER_PARSER_ELEM_TOKEN &&
        statement_elem->next->next->token.type == COMPILER_TOKEN_NAME_SYMBOL &&
        (!strcmp("[", statement_elem->next->next->token.token))) {
        sprintf(print_line, "// TODO index var\n");
        fwrite(print_line, 1, strlen(print_line), out_file);
        return 0;
    }

    if (statement_elem->next->next->next->type != COMPILER_PARSER_ELEM_EXPRESSION) {
        printf("Compiler code gen vm let statement - expected expression.\n");
        return -1;
    }

    if (compiler_code_gen_vm_expression(class_name,
                                        fun_name,
                                        routine_kind,
                                        statement_elem->next->next->next,
                                        var_class_table,
                                        var_subroutine_table,
                                        n_class_field,
                                        n_fun_var,
                                        out_file) < 0) {
        return -1;
    }

    if (statement_elem->next->type != COMPILER_PARSER_ELEM_TOKEN ||
        statement_elem->next->token.type != COMPILER_TOKEN_NAME_IDENTIFIER) {
        printf("Compiler code gen vm let statement - expected varName.\n");
        return -1;
    }
    char* var_name = statement_elem->next->token.token;

    struct compiler_code_gen_var* var = compiler_code_gen_find_fun_var(var_class_table, var_subroutine_table, var_name);
    if (var == NULL) {
        printf("Compiler code gen vm let statement - unknown varName = %s.\n", var_name);
        return -1;
    }

    if (compiler_code_gen_vm_pop_var(class_name, var, out_file) < 0)
        return -1;

    return 0;
}

int compiler_code_gen_vm_if_statement(char* class_name,
                                      char* fun_name,
                                      uint8_t routine_kind,
                                      struct compiler_parser_elem* statement_elem,
                                      struct compiler_code_var_table* var_class_table,
                                      struct compiler_code_var_table* var_subroutine_table,
                                      int n_class_field,
                                      int n_fun_var,
                                      FILE* out_file)
{
    char print_line[2048+1];

    //
    sprintf(print_line, "// if statement\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

int compiler_code_gen_vm_while_statement(char* class_name,
                                         char* fun_name,
                                         uint8_t routine_kind,
                                         struct compiler_parser_elem* statement_elem,
                                         struct compiler_code_var_table* var_class_table,
                                         struct compiler_code_var_table* var_subroutine_table,
                                         int n_class_field,
                                         int n_fun_var,
                                         FILE* out_file)
{
    char print_line[2048+1];

    //
    sprintf(print_line, "// while statement\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

int compiler_code_gen_vm_do_statement(char* class_name,
                                      char* fun_name,
                                      uint8_t routine_kind,
                                      struct compiler_parser_elem* statement_elem,
                                      struct compiler_code_var_table* var_class_table,
                                      struct compiler_code_var_table* var_subroutine_table,
                                      int n_class_field,
                                      int n_fun_var,
                                      FILE* out_file)
{
    char print_line[2048+1];

    //
    sprintf(print_line, "// do statement\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}

int compiler_code_gen_vm_return_statement(char* class_name,
                                          char* fun_name,
                                          uint8_t routine_kind,
                                          struct compiler_parser_elem* statement_elem,
                                          struct compiler_code_var_table* var_class_table,
                                          struct compiler_code_var_table* var_subroutine_table,
                                          int n_class_field,
                                          int n_fun_var,
                                          FILE* out_file)
{
    char print_line[2048+1];

    //
    sprintf(print_line, "// return statement\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    return 0;
}


int compiler_code_gen_vm_statement(char* class_name,
                                   char* fun_name,
                                   uint8_t routine_kind,
                                   struct compiler_parser_elem* statement_elem,
                                   struct compiler_code_var_table* var_class_table,
                                   struct compiler_code_var_table* var_subroutine_table,
                                   int n_class_field,
                                   int n_fun_var,
                                   FILE* out_file)
{
    // statement: letStatement | ifStatement | whileStatement | doStatement | returnStatement

    if (statement_elem == NULL) {
        printf("Compiler code gen vm subroutine - err NULL statement parse elem.\n");
        return -1;
    }

    int (*p_vm_statement) (char* class_name,
                           char* fun_name,
                           uint8_t routine_kind,
                           struct compiler_parser_elem* statement_elem,
                           struct compiler_code_var_table* var_class_table,
                           struct compiler_code_var_table* var_subroutine_table,
                           int n_class_field,
                           int n_fun_var,
                           FILE* out_file) = NULL;

    switch (statement_elem->type) {
        case COMPILER_PARSER_ELEM_LET_STATEMENT:
            p_vm_statement = compiler_code_gen_vm_let_statement;
            break;
        case COMPILER_PARSER_ELEM_IF_STATEMENT:
            p_vm_statement = compiler_code_gen_vm_if_statement;
            break;
        case COMPILER_PARSER_ELEM_WHILE_STATEMENT:
            p_vm_statement = compiler_code_gen_vm_while_statement;
            break;
        case COMPILER_PARSER_ELEM_DO_STATEMENT:
            p_vm_statement = compiler_code_gen_vm_do_statement;
            break;
        case COMPILER_PARSER_ELEM_RETURN_STATEMENT:
            p_vm_statement = compiler_code_gen_vm_return_statement;
            break;
        default:
            break;
    }

    if (p_vm_statement == NULL) {
        printf("Compiler code gen vm subroutine - err unknown %d statement type in parse elem.\n", statement_elem->type);
        return -1;
    }

    if (p_vm_statement(class_name,
                       fun_name,
                       routine_kind,
                       statement_elem,
                       var_class_table,
                       var_subroutine_table,
                       n_class_field,
                       n_fun_var,
                       out_file) < 0) {
        return -1;
    }

    return 0;
}

int compiler_code_gen_vm_statements(char* class_name,
                                    char* fun_name,
                                    uint8_t routine_kind,
                                    struct compiler_parser_elem* statements_elem,
                                    struct compiler_code_var_table* var_class_table,
                                    struct compiler_code_var_table* var_subroutine_table,
                                    int n_class_field,
                                    int n_fun_var,
                                    FILE* out_file)
{
    // statements: statement*

    char print_line[2048+1];

    //
    sprintf(print_line, "// statements\n");
    fwrite(print_line, 1, strlen(print_line), out_file);

    if (statements_elem == NULL) {
        printf("Compiler code gen vm subroutine - err NULL statements parse elem.\n");
        return -1;
    }
    if (statements_elem->type != COMPILER_PARSER_ELEM_STATEMENTS) {
        printf("Compiler code gen vm subroutine - err not statements in statement parse elem.\n");
        return -1;
    }
    statements_elem = statements_elem->child;

    while (statements_elem != NULL) {
        if (compiler_code_gen_vm_statement(class_name,
                                           fun_name,
                                           routine_kind,
                                           statements_elem,
                                           var_class_table,
                                           var_subroutine_table,
                                           n_class_field,
                                           n_fun_var,
                                           out_file) < 0) {
            return -1;
        }
        statements_elem = statements_elem->next;
    }

    return 0;
}

int compiler_code_gen_vm_routine(char* class_name,
                                 char* fun_name,
                                 uint8_t routine_kind,
                                 struct compiler_parser_elem* fun_elem,
                                 struct compiler_code_var_table* var_class_table,
                                 struct compiler_code_var_table* var_subroutine_table,
                                 int n_class_field,
                                 int n_fun_var,
                                 FILE* out_file)
{
    char print_line[2048+1];

    // function
    sprintf(print_line, "function %s.%s %d\n",
                class_name,
                fun_name,
                n_fun_var);
    fwrite(print_line, 1, strlen(print_line), out_file);

    if (COMPILER_CODE_GEN_ROUTINE_KIND_CONSTRUCTOR == routine_kind) {
        sprintf(print_line, "// constructor\n"
                            "push constant %d\n"
                            "call Memory.alloc 1\n"
                            "pop pointer 0\n",
                n_class_field);
    } else if (COMPILER_CODE_GEN_ROUTINE_KIND_METHOD == routine_kind) {
        sprintf(print_line, "// method\n"
                            "push argument 0\n"
                            "pop pointer 0\n");
    } else if (COMPILER_CODE_GEN_ROUTINE_KIND_FUNCTION == routine_kind) {
            sprintf(print_line, "// function\n");
    } else {
        sprintf(print_line, "// Unknown kind\n");
    }
    fwrite(print_line, 1, strlen(print_line), out_file);

    // subroutineBody: '{' varDec* statements '}'
    
    if (fun_elem == NULL) {
        printf("Compiler code gen vm subroutine - err NULL fun parse elem.\n");
        return -1;
    }
    if (fun_elem->type != COMPILER_PARSER_ELEM_SUBROUTINE_DEC) {
        printf("Compiler code gen vm subroutine - err not subroutineDec in fun parse elem.\n");
        return -1;
    }
    fun_elem = fun_elem->child;

    while (fun_elem != NULL) {
        if (fun_elem->type == COMPILER_PARSER_ELEM_SUBROUTINE_BODY)
               break;
        fun_elem = fun_elem->next;
    }
    if (fun_elem == NULL) {
        printf("Compiler code gen vm subroutine - err not subroutineBody in fun parse elem.\n");
        return -1;
    }
    fun_elem = fun_elem->child;

    while (fun_elem != NULL) {
        if (fun_elem->type == COMPILER_PARSER_ELEM_STATEMENTS)
               break;
        fun_elem = fun_elem->next;
    }
    if (fun_elem == NULL) {
        printf("Compiler code gen vm subroutine - err not statements in fun parse elem.\n");
        return -1;
    }

    if (compiler_code_gen_vm_statements(class_name,
                                        fun_name,
                                        routine_kind,
                                        fun_elem,
                                        var_class_table,
                                        var_subroutine_table,
                                        n_class_field,
                                        n_fun_var,
                                        out_file) < 0) {
        return -1;
    }

    sprintf(print_line, "\n");
    fwrite(print_line, 1, strlen(print_line), out_file);
    return 0;
}

int compiler_do(char* jack_file_name, uint8_t out_format, char* out_file_name)
{
    FILE* jack_file;
    FILE* out_file;
    char print_line[2048+1];
    struct compiler_parser_elem* class_elem = NULL;

    jack_file = fopen(jack_file_name, "r");
    if (jack_file == NULL) {
        printf("compiler: err to open jack file %s\n", jack_file_name);
        return -1;
    }

    out_file = fopen(out_file_name, "w");
    if (out_file == NULL) {
        printf("compiler: err to create out file %s\n", out_file_name);
        fclose(jack_file);
        return -1;
    }

    if (out_format == COMPILER_OUT_FORMAT_XML_TOKENS) {
        sprintf(print_line, "<tokens>\n");
        fwrite(print_line, 1, strlen(print_line), out_file);
    }

    tokens_t tokens;
    compiler_init_tokens(&tokens);

    char line[COMPILER_LINE_MAX_LEN+1];
    char strip_line[COMPILER_LINE_MAX_LEN+1];
    int n = 0;

    uint8_t is_in_comment = 0;

    while (fgets(line, sizeof(line), jack_file) != NULL) {

        printf("compiler: line %d <%s>, len = %d\n", n, line, (int)strlen(line));
        if (strlen(line) >= COMPILER_LINE_MAX_LEN) {
            printf("compiler: -err line %d - str of line is too big = %d \n", n, (int)strlen(line));
            break;
        }
        if (compiler_strip_line(line, &is_in_comment, strip_line) < 0) {
            printf("compiler: to strip line error\n");
            goto err;
        }
        printf("compiler: strip_line %d <%s>, len = %d\n", n, strip_line, (int)strlen(strip_line));
        if (strlen(strip_line) == 0)
            continue;
        if (out_format == COMPILER_OUT_FORMAT_STRIP) {
            fwrite(strip_line, 1, strlen(strip_line), out_file);
            fwrite("\n", 1, 1, out_file);
        }

        ++n;

        if (out_format == COMPILER_OUT_FORMAT_STRIP)
            continue;

        char* next = strip_line;
        while (strlen(next) != 0) {
            char token[COMPILER_TOKEN_MAX_LEN+1];
            char *token_name;
            if (compiler_get_token(&next, token) < 0)
                break;
            printf("compiler: token = %s, next = <%s>\n", token, next);
            uint8_t token_type = compiler_get_token_type(token);
            if (token_type == COMPILER_TOKEN_NAME_STRING_CONST)
                compiler_skip_string_quotes(token);
            if (out_format == COMPILER_OUT_FORMAT_XML_TOKENS) {
                sprintf(print_line, "    <%s> %s </%s>\n",
                        compiler_token_name[token_type],
                        (!strcmp(token, "<")) ? "&lt;" : token,
                        compiler_token_name[token_type]);
                fwrite(print_line, 1, strlen(print_line), out_file);
            }
            if (compiler_add_to_tokens(&tokens, token_type, token) < 0) {
                printf("compiler: err max number of tokens %d\n", tokens.num);
                goto err;
            }
        }
    }

    if (out_format == COMPILER_OUT_FORMAT_XML_TOKENS) {
        sprintf(print_line, "</tokens>");
        fwrite(print_line, 1, strlen(print_line), out_file);
    }

    if (out_format >= COMPILER_OUT_FORMAT_XML_PARSER) {
        int i = 0;
        if (compiler_parse_class(&tokens, &i, &class_elem) < 0) {
            printf("compiler: error parse class\n");
            ///TODO goto err;
        }
    }

    if (out_format == COMPILER_OUT_FORMAT_XML_PARSER)
        compiler_parser_print_elems(class_elem, out_file);

    if (out_format >= COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE) {
        token_t* class_name_token;
        if ((class_name_token = compiler_code_gen_get_class_name_token(class_elem)) == NULL)
            goto err;

        // Class table
        compiler_code_gen_init_var_table(&compiler_code_var_class_table);
        struct compiler_parser_elem* parse_elem = class_elem;
        int n_field = 0;
        int n_static = 0;
        int n_var;
        while ((n_var = compiler_gen_code_parse_class_var(&parse_elem, &compiler_code_var_class_table, &n_field, &n_static)) > 0)
            ;
        if (n_var < 0) {
            goto err;
        }

        if (out_format == COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE) {
            int i;
            struct compiler_code_gen_var* var;
            
            sprintf(print_line, "class %s symbol table:\n", class_name_token->token);
            fwrite(print_line, 1, strlen(print_line), out_file);
            i = 0;
            while ((var = compiler_code_gen_get_from_var_table(&compiler_code_var_class_table, i)) != NULL) {
                sprintf(print_line, "\t%3d %32s %16s %10s %8d\n", i, 
                            var->name,
                            var->type,
                            compiler_code_gen_var_ar_kind_name[var->kind],
                            var->num);
                fwrite(print_line, 1, strlen(print_line), out_file);
                ++i;
            }
        }
        
        if (out_format == COMPILER_OUT_FORMAT_VM) {
                // File, class name
                sprintf(print_line, "// compile code gen, fileName=%s, className=%s\n",
                              jack_file_name,
                              class_name_token->token);
                fwrite(print_line, 1, strlen(print_line), out_file);
                //
        }
        
        // Function tables.
        token_t* fun_name_token;
        struct compiler_parser_elem* fun_elem = class_elem;
        while ((fun_name_token = compiler_code_gen_get_next_fun(&fun_elem)) != NULL) {
            //
            uint8_t routine_kind = compiler_code_gen_get_routine_kind(fun_elem);
            //
            if (out_format == COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE) {
                sprintf(print_line, "\n\troutine %10s %16s() symbol table:\n",
                            compiler_code_gen_routine_ar_kind_name[routine_kind],
                            fun_name_token->token);
                fwrite(print_line, 1, strlen(print_line), out_file);
            }

            compiler_code_gen_init_var_table(&compiler_code_var_subroutine_table);
            // Arg vars
            int n_arg;
            if ((n_arg = compiler_gen_code_parse_fun_arg_var(fun_elem, &compiler_code_var_subroutine_table)) < 0) {
                goto err;
            }
            // Local vars
            parse_elem = fun_elem;
            int n_local = 0;
            while ((n_var = compiler_gen_code_parse_fun_var(&parse_elem, &compiler_code_var_subroutine_table, &n_local)) > 0)
                ;
            if (n_local < 0) {
                goto err;
            }

            if (out_format == COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE) {
                int i = 0;
                struct compiler_code_gen_var* var;
                while ((var = compiler_code_gen_get_from_var_table(&compiler_code_var_subroutine_table, i)) != NULL) {
                    sprintf(print_line, "\t\t%3d %32s %16s %10s %8d\n", i,
                              var->name,
                              var->type,
                              compiler_code_gen_var_ar_kind_name[var->kind],
                              var->num);
                     fwrite(print_line, 1, strlen(print_line), out_file);
                     ++i;
                }
            }

            if (out_format == COMPILER_OUT_FORMAT_VM) {
                // Routine
                if (compiler_code_gen_vm_routine(class_name_token->token,
                                                 fun_name_token->token,
                                                 routine_kind,
                                                 fun_elem,
                                                 &compiler_code_var_class_table,
                                                 &compiler_code_var_subroutine_table,
                                                 n_field,
                                                 n_local,
                                                 out_file) < 0) {
                    goto err;
                }
            }
        }
    }


err:
    fclose(jack_file);
    fclose(out_file);
    return 0;
}


int main()
{
	printf("Test:\n");
    computer_t computer;


	dff_t dff; 
        dff16_t dff16;

	dff_init(&dff);
	dff16_init(&dff16);

	bit_t in, out;
	bit16_t in_bit16, out_bit16;

	in = 1;
	in_bit16 = bit16_ones;

	print_bit("Set in", in);
	print_bit16("Set in16", &in_bit16);
	out = dff_get_output(&dff);
	out_bit16 = dff16_get_output(&dff16);
	print_dff(&dff);
	print_dff16(&dff16);
	printf("Run\n");
	dff_run(&dff, in);
	dff16_run(&dff16, &in_bit16);
	out = dff_get_output(&dff);
	out_bit16 = dff16_get_output(&dff16);
	print_dff(&dff);
	print_dff16(&dff16);


        bit_t a,b;
        a = 0; b = 0; out = 0;
        nand(a, b, &out);
        printf("a = %s, b = %s, nand_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 0; out = 0;
        nand(a, b, &out);
        printf("a = %s, b = %s, nand_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 0; b = 1; out = 0;
        nand(a, b, &out);
        printf("a = %s, b = %s, nand_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 1; out = 0;
        nand(a, b, &out);
        printf("a = %s, b = %s, nand_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        
	a = 0; b = 0; out = 0;
        gate_and(a, b, &out);
        printf("a = %s, b = %s, and_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 0; out = 0;
        gate_and(a, b, &out);
        printf("a = %s, b = %s, and_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 0; b = 1; out = 0;
        gate_and(a, b, &out);
        printf("a = %s, b = %s, and_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 1; out = 0;
        gate_and(a, b, &out);
        printf("a = %s, b = %s, and_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        
	a = 0; b = 0; out = 0;
	gate_or(a, b, &out);
        printf("a = %s, b = %s, or_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 0; out = 0;
        gate_or(a, b, &out);
        printf("a = %s, b = %s, or_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 0; b = 1; out = 0;
        gate_or(a, b, &out);
        printf("a = %s, b = %s, or_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 1; out = 0;
        gate_or(a, b, &out);
        printf("a = %s, b = %s, or_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 

	a = 0; b = 0; out = 0;
	num_not = 0; num_or = 0; num_and = 0;
        gate_xor(a, b, &out);
	printf("Xor: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
        printf("a = %s, b = %s, xor_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 0; out = 0;
        gate_xor(a, b, &out);
        printf("a = %s, b = %s, xor_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 0; b = 1; out = 0;
        gate_xor(a, b, &out);
        printf("a = %s, b = %s, xor_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
        a = 1; b = 1; out = 0;
        gate_xor(a, b, &out);
        printf("a = %s, b = %s, xor_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(out)); 
	
	
        bit_t sel = 0;	
	a = 0; b = 0; out = 0;
	num_not = 0; num_or = 0; num_and = 0;
        mux(a, b, sel, &out);
	printf("Mux: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 1; b = 0; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 0; b = 1; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 1; b = 1; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
        sel = 1;	
	a = 0; b = 0; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 1; b = 0; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 0; b = 1; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
	a = 1; b = 1; out = 0;
        mux(a, b, sel, &out);
        printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sel), bit_to_string(out)); 
       
        {
		bit16_t a = bit16_zeros;
		bit16_t b = bit16_ones;
		bit_t sel = 0;
		bit16_t out = bit16_zeros;
		char str1[BIT16_SIZE];
		char str2[BIT16_SIZE];
		char str3[BIT16_SIZE];

	        num_not = 0; num_or = 0; num_and = 0;
		mux16(&a, &b, sel, &out);
	        printf("Mux16: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
                printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit16_to_string(str1, &a), bit16_to_string(str2, &b), bit_to_string(sel), bit16_to_string(str3, &out)); 
		
		sel = 1;
		mux16(&a, &b, sel, &out);
                printf("a = %s, b = %s, sel - %s, mux_out = %s\n", bit16_to_string(str1, &a), bit16_to_string(str2, &b), bit_to_string(sel), bit16_to_string(str3, &out)); 

	}

	{
		bit_t a, b, c, carry, sum;

		a = 0; b = 0; carry = 0; sum = 0;
	        num_not = 0; num_or = 0; num_and = 0;
		hadd(a, b, &sum, &carry);
	        printf("hadd: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
		printf("a = %s, b = %s, haddr.sum = %s, haddr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 0; carry = 0; sum = 0;
		hadd(a, b, &sum, &carry);
		printf("a = %s, b = %s, haddr.sum = %s, haddr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sum), bit_to_string(carry));
		a = 0; b = 1; carry = 0; sum = 0;
		hadd(a, b, &sum, &carry);
		printf("a = %s, b = %s, haddr.sum = %s, haddr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 1; carry = 0; sum = 0;
		hadd(a, b, &sum, &carry);
		printf("a = %s, b = %s, haddr.sum = %s, haddr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(sum), bit_to_string(carry));

		a = 0; b = 0, c = 0; carry = 0; sum = 0;
	        num_not = 0; num_or = 0; num_and = 0;
		add(a, b, c, &sum, &carry);
	        printf("add: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 0, c = 0; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 0; b = 1, c = 0; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 1, c = 0; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 0; b = 0, c = 1; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 0, c = 1; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 0; b = 1, c = 1; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
		a = 1; b = 1, c = 1; carry = 0; sum = 0;
		add(a, b, c, &sum, &carry);
		printf("a = %s, b = %s, c = %s, addr.sum = %s, addr.carry = %s\n", bit_to_string(a), bit_to_string(b), bit_to_string(c), bit_to_string(sum), bit_to_string(carry));
	
	
	
	}

	{
		bit16_t bit16 = bit16_ones;
		int16_t val = 0, val2 = 0;
		char str[BIT16_SIZE];


		val = 1;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = -1;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = 101;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = 16536;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = 16535;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = 4095;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		val = 4096;
		int16_to_bit16(val, &bit16);
		printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));


		{		
			bit16_t bit16 = {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; 
			val = bit16_to_int16(&bit16);         
			printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		}
		{		
			bit16_t bit16 = {{1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}; 
			val = bit16_to_int16(&bit16);         
			printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		}
		{		
			bit16_t bit16 = {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}; 
			val = bit16_to_int16(&bit16);         
			printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		}
		{		
			bit16_t bit16 = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}; 
			val = bit16_to_int16(&bit16);         
			printf("val = %d, bit16 = %s\n", val, bit16_to_string(str, &bit16));
		}

	}
	
	{
		bit16_t in16;
		bit_t load;
		bit16_t out16;
		reg16_t reg16;
	        num_dff = 0; num_not = 0; num_or = 0; num_and = 0;
		reg16_init(&reg16, "A16");


		char str1[BIT16_SIZE];
		char str2[BIT16_SIZE];

		int16_to_bit16(1, &in16);
		load = 0;
		out16 = reg16_get_output(&reg16);
		printf("in16 = %s, load = %s, reg16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		reg16_run(&reg16, &in16, load);
	        printf("reg16: num not=%u, or=%u, and=%u, dff=%u\n", num_not, num_or, num_and, num_dff);
		out16 = reg16_get_output(&reg16);
		printf("in16 = %s, load = %s, reg16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(2, &in16);
		load = 1;
		out16 = reg16_get_output(&reg16);
		printf("in16 = %s, load = %s, reg16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		reg16_run(&reg16, &in16, load);
		out16 = reg16_get_output(&reg16);
		printf("in16 = %s, load = %s, reg16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(load), bit16_to_string(str2, &out16));
		

	}

	{
		bit16_t in16;
		bit_t reset, load;
		bit16_t out16;
		pc16_t pc16;
	        num_dff = 0; num_not = 0; num_or = 0; num_and = 0;
		pc16_init(&pc16);


		char str1[BIT16_SIZE];
		char str2[BIT16_SIZE];

		int16_to_bit16(1, &in16);
		reset = 1; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
	        printf("pc16: num not=%u, or=%u, and=%u, dff=%u\n", num_not, num_or, num_and, num_dff);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(1, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(1, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(1, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 1;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 1; load = 1;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		
		int16_to_bit16(55, &in16);
		reset = 0; load = 0;
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		printf("Run.\n");
		pc16_run(&pc16, reset, &in16, load);
		out16 = pc16_get_output(&pc16);
		printf("in16 = %s, reset = %s, load = %s, pc16.out = %s\n", bit16_to_string(str1, &in16), bit_to_string(reset), bit_to_string(load), bit16_to_string(str2, &out16));
		print_pc16(&pc16);

	}

	{
		bit16_t x16, y16;
		bit_t zx, nx, zy, ny, f, nf, nz, ng;
		bit16_t out16;
	        num_not = 0; num_or = 0; num_and = 0;

		char str1[BIT16_SIZE];
		char str2[BIT16_SIZE];
		char str3[BIT16_SIZE];

		int16_to_bit16(2222, &x16);
		int16_to_bit16(5555, &y16);
		zx = 0; nx = 0; zy = 0; ny = 0; f = 0; nf = 0;

		alu(&x16, &y16, zx, nx, zy, ny, f, nf, &out16, &nz, &ng);	
		printf("alu: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
		printf("x16 = %s (%d), y16 = %s (%d), zx = %s, nx = %s, zy = %s, ny = %s, f = %s, nf = %s, alu.out = %s (%d), alu.nz = %s, alu.ng = %s\n",
				bit16_to_string(str1, &x16), bit16_to_int16(&x16), bit16_to_string(str2, &y16), bit16_to_int16(&y16), bit_to_string(zx), bit_to_string(nx),
				bit_to_string(zy), bit_to_string(ny), bit_to_string(f), bit_to_string(nf),
				bit16_to_string(str3, &out16), bit16_to_int16(&out16), bit_to_string(nz), bit_to_string(ng));
	
	        int16_to_bit16(3333, &x16);
                int16_to_bit16(8888, &y16);
                zx = 0; nx = 0; zy = 1; ny = 1; f = 0; nf = 0;

                alu(&x16, &y16, zx, nx, zy, ny, f, nf, &out16, &nz, &ng);
                printf("alu: num not=%u, or=%u, and=%u\n", num_not, num_or, num_and);
                printf("x16 = %s (%d), y16 = %s (%d), zx = %s, nx = %s, zy = %s, ny = %s, f = %s, nf = %s, alu.out = %s (%d), alu.nz = %s, alu.ng = %s\n",
                                bit16_to_string(str1, &x16), bit16_to_int16(&x16), bit16_to_string(str2, &y16), bit16_to_int16(&y16), bit_to_string(zx), bit_to_string(nx),
                                bit_to_string(zy), bit_to_string(ny), bit_to_string(f), bit_to_string(nf),
                                bit16_to_string(str3, &out16), bit16_to_int16(&out16), bit_to_string(nz), bit_to_string(ng));
	}
	
	
	{
		bit_t reset = 0;
		
	    num_dff = 0; num_not = 0; num_or = 0; num_and = 0;

		char str1[BIT16_SIZE];
		char str2[BIT16_SIZE];
		char str3[BIT16_SIZE];
		char str4[BIT16_SIZE];
		char str5[BIT16_SIZE];
		char str_instruction[INSTRUCTION_SIZE];

	        memory_debug = 1;
	        screen_debug = 1;

	    computer_init(&computer);
	    
#if 1
	    compiler_do("Main.jack", COMPILER_OUT_FORMAT_STRIP, "Main.jack_strip");
	    compiler_do("Main.jack", COMPILER_OUT_FORMAT_XML_TOKENS, "MainT.xml");
	    compiler_do("Main.jack", COMPILER_OUT_FORMAT_XML_PARSER, "Main.xml");
	    compiler_do("Main.jack", COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE, "Main.symbol_table");
	    compiler_do("Main.jack", COMPILER_OUT_FORMAT_VM, "Main.vm");
#endif
	    compiler_do("Square.jack", COMPILER_OUT_FORMAT_STRIP, "Square.jack_strip");
	    compiler_do("Square.jack", COMPILER_OUT_FORMAT_XML_TOKENS, "SquareT.xml");
	    compiler_do("Square.jack", COMPILER_OUT_FORMAT_XML_PARSER, "Square.xml");
	    compiler_do("Square.jack", COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE, "Square.symbol_table");
	    compiler_do("Square.jack", COMPILER_OUT_FORMAT_VM, "Square.vm");
#if 1
	    compiler_do("SquareGame.jack", COMPILER_OUT_FORMAT_STRIP, "SquareGame.jack_strip");
	    compiler_do("SquareGame.jack", COMPILER_OUT_FORMAT_XML_TOKENS, "SquareGameT.xml");
	    compiler_do("SquareGame.jack", COMPILER_OUT_FORMAT_XML_PARSER, "SquareGame.xml");
	    compiler_do("SquareGame.jack", COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE, "SquareGame.symbol_table");
	    compiler_do("SquareGame.jack", COMPILER_OUT_FORMAT_VM, "SquareGame.vm");
	    
	    compiler_do("Main2.jack", COMPILER_OUT_FORMAT_STRIP, "Main2.jack_strip");
	    compiler_do("Main2.jack", COMPILER_OUT_FORMAT_XML_TOKENS, "Main2T.xml");
	    compiler_do("Main2.jack", COMPILER_OUT_FORMAT_XML_PARSER, "Main2.xml");
	    compiler_do("Main2.jack", COMPILER_OUT_FORMAT_VAR_SYMBOL_TABLE, "Main2.symbol_table");
	    compiler_do("Main2.jack", COMPILER_OUT_FORMAT_VM, "Main2.vm");
#endif
	    // TODO
	    return 0;
	
	    vm_translate("fibonacci.vm", VM_OUT_FORMAT_STRIP, "fibonacci.vm_strip");
	    vm_translate("fibonacci.vm", VM_OUT_FORMAT_ASM, "fibonacci.asm");
	    
	    asm_symbol_tbl_t symbol_tbl;
	    asm_init_symbol_tbl(&symbol_tbl);
	    asm_translate("fibonacci.asm", ASM_OUT_FORMAT_SYMBOLS, &symbol_tbl, "fibonacci.asm_symbol");
	    asm_translate("fibonacci.asm", ASM_OUT_FORMAT_LITERALS, &symbol_tbl, "fibonacci.asm_literal");
	    asm_translate("fibonacci.asm", ASM_OUT_FORMAT_BIN, &symbol_tbl, "fibonacci.hack");

        int res = rom32k_load(&(computer.rom32k), "fibonacci.hack");
        if (res != 0) {
            return -1;
        }

        // R0 = 5;
        computer.memory.ram16k.mem[0]=15;

#if 0
		int i;
		for (i = 0; i < 500/*AR_SIZEOF(code)*/; i++)
        {
        pc = pc16_get_output(&(computer.cpu.pc16));
        instruction = rom32k_get_output(&(computer.rom32k), &pc);
	    in_m = memory_get_output(&(computer.memory));
		printf("Run.\n");
		
		print_reg16(&(computer.cpu.reg_a16));
		print_reg16(&(computer.cpu.reg_d16));
		printf("in_m = %s (0x%x, %d), instruction = %s (0x%x, %d), reset = %s,\n"
		       "instruction: %s\n"
		        "cpu.out_m = %s (0x%x, %d), cpu.wr_m =%s, cpu.addr_m = %s (0x%x), cpu.pc = %s (0x%x)\n",
		        bit16_to_string(str1, &in_m), bit16_to_int16(&in_m), bit16_to_int16(&in_m), 
		        bit16_to_string(str2, &instruction), bit16_to_int16(&instruction), bit16_to_int16(&instruction), bit_to_string(reset),
		        instruction_bit16_to_asm(str_instruction, &instruction),
		        bit16_to_string(str3, &out_m), bit16_to_int16(&out_m), bit16_to_int16(&out_m),
		        bit_to_string(wr_m), bit16_to_string(str4, &addr_m), bit16_to_int16(&addr_m), bit16_to_string(str5, &pc), bit16_to_int16(&pc));
	    }
	    ram16k_dump(&(computer.memory.ram16k), 0, 10);
	    screen_dump(&(computer.memory.screen), 0, 40, 0, 80);
#endif
	}
	
	
	const int screenWidth = 512+512+64;
    const int screenHeight = 256+256;

    SDL_Init(SDL_INIT_EVERYTHING);
    SDL_Window* window = SDL_CreateWindow("HAck", SDL_WINDOWPOS_UNDEFINED, 
                                                  SDL_WINDOWPOS_UNDEFINED, screenWidth, screenHeight, SDL_WINDOW_SHOWN);
    
    SDL_Renderer* render = SDL_CreateRenderer(window, -1, 0);
    SDL_SetRenderDrawColor(render, 220, 220, 220, 255);
    //SDL_UpdateWindowSurface(window);
    TTF_Init();
    TTF_Font* Sans24 = TTF_OpenFont("TTF/FreeMono.ttf", 24);
    TTF_Font* Sans18 = TTF_OpenFont("TTF/FreeMono.ttf", 18);
    
    SDL_Rect rect;
    rect.x = 512/2;
    rect.y = 256/2;
    rect.w = 512;
    rect.h = 256;
    SDL_Color color = {0, 0, 255};
    SDL_Color color_ram   = {255, 0, 0};
    SDL_Color color_kbd   = {0, 0, 255};
    SDL_Color color_stack = {0, 0, 255};
    
    bit16_t bit16;
    bit16_t bit16_instr;
    int16_t bit16_num;
    char str[256];
    char asm_str[100];
    char str_ram[1024];

    SDL_Surface* image;
    SDL_Texture* texture1;

    SDL_Surface* surface_pc16; 
    SDL_Texture* message_pc16;
    SDL_Surface* surface_a16; 
    SDL_Texture* message_a16;
    SDL_Surface* surface_d16; 
    SDL_Texture* message_d16;
    SDL_Surface* surface_kbd; 
    SDL_Texture* message_kbd;
    
    SDL_Rect reg_rect;
    
    
    SDL_Surface* surface_ram;
    SDL_Texture* message_ram;
    SDL_Surface* surface_ram_reg;
    SDL_Texture* message_ram_reg;
    SDL_Surface* surface_ram_stack;
    SDL_Texture* message_ram_stack;
    
    
    SDL_Event event;
    const Uint8* m_keystates;

    uint8_t is_run = 0;

    while (1) {
    
        uint8_t is_break = 0;
    
        bit16 = pc16_get_output(&(computer.cpu.pc16));
        bit16_instr = rom32k_get_output(&(computer.rom32k), &bit16);
        instruction_bit16_to_asm(asm_str, &bit16_instr);
        bit16_num = bit16_to_int16(&bit16);
        sprintf(str, "PC: 0x%.4x    %s", bit16_num, asm_str);
        surface_pc16 = TTF_RenderText_Solid(Sans24, str, color); 
        message_pc16 = SDL_CreateTextureFromSurface(render, surface_pc16);

        bit16 = reg16_get_output(&(computer.cpu.reg_a16));
        bit16_num = bit16_to_int16(&bit16);
        sprintf(str, " A: 0x%.4x", bit16_num);
        surface_a16 = TTF_RenderText_Solid(Sans24, str, color); 
        message_a16 = SDL_CreateTextureFromSurface(render, surface_a16);
        
        bit16 = reg16_get_output(&(computer.cpu.reg_d16));
        bit16_num = bit16_to_int16(&bit16);
        sprintf(str, " D: 0x%.4x", bit16_num);
        surface_d16 = TTF_RenderText_Solid(Sans24, str, color); 
        message_d16 = SDL_CreateTextureFromSurface(render, surface_d16);

        strcpy(str_ram, "RAM:\n");
        int i;
        for (i = 16; i < 16+48; i+=4) {
            sprintf(str, "%.2x: %.4x %.4x %.4x %.4x\n", i, 
                          (uint16_t)computer.memory.ram16k.mem[i],
                          (uint16_t)computer.memory.ram16k.mem[i+1],
                          (uint16_t)computer.memory.ram16k.mem[i+2],
                          (uint16_t)computer.memory.ram16k.mem[i+3]);
            strcat(str_ram, str);
        }
        surface_ram = TTF_RenderText_Solid_Wrapped(Sans18, str_ram, color_ram, 256+64);
        message_ram = SDL_CreateTextureFromSurface(render, surface_ram);

        strcpy(str_ram, "STACK: +256 (0x100)\n");
        for (i = 0; i < 48; i+=4) {
            sprintf(str, "%.2x: %.4x %.4x %.4x %.4x\n", i, 
                          (uint16_t)computer.memory.ram16k.mem[256+i],
                          (uint16_t)computer.memory.ram16k.mem[256+i+1],
                          (uint16_t)computer.memory.ram16k.mem[256+i+2],
                          (uint16_t)computer.memory.ram16k.mem[256+i+3]);
            strcat(str_ram, str);
        }
        surface_ram_stack = TTF_RenderText_Solid_Wrapped(Sans18, str_ram, color_stack, 256+64);
        message_ram_stack = SDL_CreateTextureFromSurface(render, surface_ram_stack);

        sprintf(str_ram, " SP : 0x%.4x\n"
                         " LCL: 0x%.4x\n"
                         " ARG: 0x%.4x\n"
                         "THIS: 0x%.4x\n"
                         "THAT: 0x%.4x\n", 
                (uint16_t)computer.memory.ram16k.mem[0],
                (uint16_t)computer.memory.ram16k.mem[1],
                (uint16_t)computer.memory.ram16k.mem[2],
                (uint16_t)computer.memory.ram16k.mem[3],
                (uint16_t)computer.memory.ram16k.mem[4]);
        for (i = 5; i < 16; i++) {
            sprintf(str, " R%2d: 0x%.4x\n", i, (uint16_t)computer.memory.ram16k.mem[i]);
            strcat(str_ram, str);
        }
        surface_ram_reg = TTF_RenderText_Solid_Wrapped(Sans18, str_ram, color_ram, 128+64);
        message_ram_reg = SDL_CreateTextureFromSurface(render, surface_ram_reg);

        bit16 = keyboard_get_output(&(computer.memory.keyboard));
        bit16_num = bit16_to_int16(&bit16);
        sprintf(str, " KBD: 0x%.2x", bit16_num);
        surface_kbd = TTF_RenderText_Solid(Sans24, str, color_kbd); 
        message_kbd = SDL_CreateTextureFromSurface(render, surface_kbd);

        screen_save_bmp(&(computer.memory.screen), "screen.bmp");
        image = SDL_LoadBMP("screen.bmp");
        texture1 = SDL_CreateTextureFromSurface(render, image);


    
        SDL_RenderClear(render);
        SDL_RenderCopy(render, texture1, NULL, &rect);
        reg_rect.x = 8; reg_rect.y = 0; reg_rect.w = surface_pc16->w; reg_rect.h = surface_pc16->h;
        SDL_RenderCopy(render, message_pc16, NULL, &reg_rect);
        reg_rect.x = 8; reg_rect.y = surface_pc16->h; reg_rect.w = surface_a16->w; reg_rect.h = surface_a16->h;
        SDL_RenderCopy(render, message_a16, NULL, &reg_rect);
        reg_rect.x = 8; reg_rect.y = surface_pc16->h + surface_a16->h; reg_rect.w = surface_d16->w; reg_rect.h = surface_d16->h;
        SDL_RenderCopy(render, message_d16, NULL, &reg_rect);
        reg_rect.x = 64; reg_rect.y = 128-16; reg_rect.w = surface_ram_reg->w; reg_rect.h = surface_ram_reg->h;
        SDL_RenderCopy(render, message_ram_reg, NULL, &reg_rect);
        reg_rect.x = 256+512+32+16+4; reg_rect.y = surface_pc16->h; reg_rect.w = surface_ram->w; reg_rect.h = surface_ram->h;
        SDL_RenderCopy(render, message_ram, NULL, &reg_rect);
        reg_rect.x = 256+512+32+16+4;  reg_rect.y = surface_pc16->h + surface_ram->h + 8; reg_rect.w = surface_ram_stack->w; reg_rect.h = surface_ram_stack->h;
        SDL_RenderCopy(render, message_ram_stack, NULL, &reg_rect);
        reg_rect.x = 256+256-surface_kbd->w/2; reg_rect.y = 128+256+8; reg_rect.w = surface_kbd->w; reg_rect.h = surface_kbd->h;
        SDL_RenderCopy(render, message_kbd, NULL, &reg_rect);
        
        
        SDL_RenderPresent(render);

        SDL_Delay(1);
        
        while (SDL_PollEvent(&event)) {

            if (event.type == SDL_QUIT)
               is_break = 1;
            
            if (event.type == SDL_KEYDOWN) {
                
                if (event.key.keysym.scancode == SDL_SCANCODE_SPACE) {
                    computer_run(&computer, 0);
                }
                if (event.key.keysym.scancode == SDL_SCANCODE_ESCAPE) {
                    memory_init(&(computer.memory));
                    computer_run(&computer, 1);
                }
                if (event.key.keysym.scancode == SDL_SCANCODE_RETURN) {
                    is_run = !is_run;
                }
                
                keyboard_set_code(&(computer.memory.keyboard), (uint16_t)event.key.keysym.scancode);
            }
            
            if (event.type == SDL_KEYUP) {
                keyboard_set_code(&(computer.memory.keyboard), 0);
            }
        }

        if (is_run) {
            computer_run(&computer, 0);
        }
        
        SDL_DestroyTexture(texture1);
        SDL_FreeSurface(image);
        SDL_DestroyTexture(message_pc16);
        SDL_FreeSurface(surface_pc16);
        SDL_DestroyTexture(message_a16);
        SDL_FreeSurface(surface_a16);
        SDL_DestroyTexture(message_d16);
        SDL_FreeSurface(surface_d16);
        SDL_DestroyTexture(message_ram);
        SDL_FreeSurface(surface_ram);
        SDL_DestroyTexture(message_kbd);
        SDL_FreeSurface(surface_kbd);
        SDL_DestroyTexture(message_ram_reg);
        SDL_FreeSurface(surface_ram_reg);
        SDL_DestroyTexture(message_ram_stack);
        SDL_FreeSurface(surface_ram_stack);
        
        if (is_break)
            break;
    }

    SDL_DestroyRenderer(render);
    SDL_DestroyWindow(window);
    SDL_Quit();
	
	printf("==End.\n");
	return 0;
}

