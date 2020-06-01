/**
 * Macros to manipulate boundary conditions bitmask
 * @file BCMacros3D.h
 * @author Maciej Kubat (m.j.kubat@cranfield.ac.uk) 
 * @details The bitmask looks the following way
 * @verbatim
   | 70 | 68 | 66 | 64 | 62 | 59 |  58  |   57   |  54  |  51    |   48   |   45   |   42   |   39   |    36   |   33  |    30  |   27   |   24   |   21   |   18   |   15   |   12   |    9   |    6   |    3   |    0   |
   |-------------------|---------|------|--------|------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
   | boundary ID (8bit)| RESERVE | FLUID| CORNER | MAIN | DIR:18 | DIR:17 | DIR:16 | DIR:15 | DIR:14 | DIR:13 | DIR:12 | DIR:11 | DIR:10 | DIR:9  |  DIR:8 |  DIR:7 |  DIR:6 |  DIR:5 |  DIR:4 |  DIR:3 |  DIR:2 |  DIR:1 | @endverbatim
 */
#ifndef BC_MACROS3D_H
#define BC_MACROS3D_H



#define  BC3D_WALL_1   0X000000000000001
#define  BC3D_INLT_1   0X000000000000002
#define  BC3D_OUTL_1   0X000000000000003
#define  BC3D_CYCL_1   0X000000000000004
#define  BC3D_SYMP_1   0X000000000000005

#define  BC3D_WALL_2   0X000000000000008
#define  BC3D_INLT_2   0X000000000000010
#define  BC3D_OUTL_2   0X000000000000018
#define  BC3D_CYCL_2   0X000000000000020
#define  BC3D_SYMP_2   0X000000000000028

#define  BC3D_WALL_3   0X000000000000040
#define  BC3D_INLT_3   0X000000000000080
#define  BC3D_OUTL_3   0X0000000000000C0
#define  BC3D_CYCL_3   0X000000000000100
#define  BC3D_SYMP_3   0X000000000000140

#define  BC3D_WALL_4   0X000000000000200
#define  BC3D_INLT_4   0X000000000000400
#define  BC3D_OUTL_4   0X000000000000600
#define  BC3D_CYCL_4   0X000000000000800
#define  BC3D_SYMP_4   0X000000000000A00

#define  BC3D_WALL_5   0X000000000001000
#define  BC3D_INLT_5   0X000000000002000
#define  BC3D_OUTL_5   0X000000000003000
#define  BC3D_CYCL_5   0X000000000004000
#define  BC3D_SYMP_5   0X000000000005000

#define  BC3D_WALL_6   0X000000000008000
#define  BC3D_INLT_6   0X000000000010000
#define  BC3D_OUTL_6   0X000000000018000
#define  BC3D_CYCL_6   0X000000000020000
#define  BC3D_SYMP_6   0X000000000028000

#define  BC3D_WALL_7   0X000000000040000
#define  BC3D_INLT_7   0X000000000080000
#define  BC3D_OUTL_7   0X0000000000C0000
#define  BC3D_CYCL_7   0X000000000100000
#define  BC3D_SYMP_7   0X000000000140000

#define  BC3D_WALL_8   0X000000000200000
#define  BC3D_INLT_8   0X000000000400000
#define  BC3D_OUTL_8   0X000000000600000
#define  BC3D_CYCL_8   0X000000000800000
#define  BC3D_SYMP_8   0X000000000A00000

#define  BC3D_WALL_9   0X000000001000000
#define  BC3D_INLT_9   0X000000002000000
#define  BC3D_OUTL_9   0X000000003000000
#define  BC3D_CYCL_9   0X000000004000000
#define  BC3D_SYMP_9   0X000000005000000

#define  BC3D_WALL_10   0X000000008000000
#define  BC3D_INLT_10   0X000000010000000
#define  BC3D_OUTL_10   0X000000018000000
#define  BC3D_CYCL_10   0X000000020000000
#define  BC3D_SYMP_10   0X000000028000000

#define  BC3D_WALL_11   0X000000040000000
#define  BC3D_INLT_11   0X000000080000000
#define  BC3D_OUTL_11   0X0000000C0000000
#define  BC3D_CYCL_11   0X000000100000000
#define  BC3D_SYMP_11   0X000000140000000

#define  BC3D_WALL_12   0X000000200000000
#define  BC3D_INLT_12   0X000000400000000
#define  BC3D_OUTL_12   0X000000600000000
#define  BC3D_CYCL_12   0X000000800000000
#define  BC3D_SYMP_12   0X000000A00000000

#define  BC3D_WALL_13   0X000001000000000
#define  BC3D_INLT_13   0X000002000000000
#define  BC3D_OUTL_13   0X000003000000000
#define  BC3D_CYCL_13   0X000004000000000
#define  BC3D_SYMP_13   0X000005000000000

#define  BC3D_WALL_14   0X000008000000000
#define  BC3D_INLT_14   0X000010000000000
#define  BC3D_OUTL_14   0X000018000000000
#define  BC3D_CYCL_14   0X000020000000000
#define  BC3D_SYMP_14   0X000028000000000

#define  BC3D_WALL_15   0X000040000000000
#define  BC3D_INLT_15   0X000080000000000
#define  BC3D_OUTL_15   0X0000C0000000000
#define  BC3D_CYCL_15   0X000100000000000
#define  BC3D_SYMP_15   0X000140000000000

#define  BC3D_WALL_16   0X000200000000000
#define  BC3D_INLT_16   0X000400000000000
#define  BC3D_OUTL_16   0X000600000000000
#define  BC3D_CYCL_16   0X000800000000000
#define  BC3D_SYMP_16   0X000A00000000000

#define  BC3D_WALL_17   0X001000000000000
#define  BC3D_INLT_17   0X002000000000000
#define  BC3D_OUTL_17   0X003000000000000
#define  BC3D_CYCL_17   0X004000000000000
#define  BC3D_SYMP_17   0X005000000000000

#define  BC3D_WALL_18   0X008000000000000
#define  BC3D_INLT_18   0X010000000000000
#define  BC3D_OUTL_18   0X018000000000000
#define  BC3D_CYCL_18   0X020000000000000
#define  BC3D_SYMP_18   0X028000000000000

#define BC3D_WALL_B     0X040000000000000 ///< BC3D: main (001) wall
#define BC3D_INLT_B     0x080000000000000 ///< BC3D: main (010) inlet
#define BC3D_OUTL_B     0x0C0000000000000 ///< BC3D: main (011) oulet
#define BC3D_CYCL_B     0x100000000000000 ///< BC3D: main (100) periodic
#define BC3D_SYMP_B     0x140000000000000 ///< BC3D: main (101) symmetry
#define BC3D_ALL_B      0x1C0000000000000 ///< BC3D: main (111) ALL

#define BC3D_CORNER     0x200000000000000 ///< BC3D: corner
#define BC3D_FLUID      0x400000000000000 ///< BC3D: fluid

#define BC3D_NONE     0x0 ///< BC3D: type none    	(000)
#define BC3D_WALL     0x1 ///< BC3D: type wall    	(001)
#define BC3D_INLT     0x2 ///< BC3D: type inlet   	(010)
#define BC3D_OUTL     0x3 ///< BC3D: type outlet   	(011)
#define BC3D_CYCL     0x4 ///< BC3D: type periodic	(100)
#define BC3D_SYMP     0x5 ///< BC3D: type symmetry  (101)
#define BC3D_ALL      0x7 ///< BC3D: type all/any 	(111)

#define	 BC3D_1(i)	 ((i)<<0)	///<	BC3D:	set	type	to	east			@param	i	BC3D	type
#define	 BC3D_2(i)	 ((i)<<3)	///<	BC3D:	set	type	to	west			@param	i	BC3D	type
#define	 BC3D_3(i)	 ((i)<<6)	///<	BC3D:	set	type	to	north			@param	i	BC3D	type
#define	 BC3D_4(i)	 ((i)<<9)	///<	BC3D:	set	type	to	south			@param	i	BC3D	type
#define	 BC3D_5(i)	((i)<<12)	///<	BC3D:	set	type	to	top				@param	i	BC3D	type
#define	 BC3D_6(i)	((i)<<15)	///<	BC3D:	set	type	to	bottom			@param	i	BC3D	type
#define	 BC3D_7(i)	((i)<<18)	///<	BC3D:	set	type	to	north-east		@param	i	BC3D	type
#define	 BC3D_8(i)	((i)<<21)	///<	BC3D:	set	type	to	north-west		@param	i	BC3D	type
#define	 BC3D_9(i)	((i)<<24)	///<	BC3D:	set	type	to	south-east		@param	i	BC3D	type
#define	BC3D_10(i)	((i)<<27)	///<	BC3D:	set	type	to	south-west		@param	i	BC3D	type
#define	BC3D_11(i)	((i)<<30)	///<	BC3D:	set	type	to	top-east		@param	i	BC3D	type
#define	BC3D_12(i)	((i)<<33)	///<	BC3D:	set	type	to	top-west		@param	i	BC3D	type
#define	BC3D_13(i)	((i)<<36)	///<	BC3D:	set	type	to	bottom-east		@param	i	BC3D	type
#define	BC3D_14(i)	((i)<<39)	///<	BC3D:	set	type	to	bottom-west		@param	i	BC3D	type
#define	BC3D_15(i)	((i)<<42)	///<	BC3D:	set	type	to	top-north		@param	i	BC3D	type
#define	BC3D_16(i)	((i)<<45)	///<	BC3D:	set	type	to	top-south		@param	i	BC3D	type
#define	BC3D_17(i)	((i)<<48)	///<	BC3D:	set	type	to	bottom-north	@param	i	BC3D	type
#define	BC3D_18(i)	((i)<<51)	///<	BC3D:	set	type	to	bottom-south	@param	i	BC3D	type
#define	BC3D_B(i)	((i)<<54)	///<	BC3D:	set	type	to	main			@param	i	BC3D  type
#define	BC3D_C(i)	((i)<<57)	///<	BC3D:	set	type	to	corner			@param	i	true/false
#define	BC3D_F(i)	((i)<<58)	///<	BC3D:	set	type	to	fluid			@param	i	true/false
/**
 * @brief BC3D: set type to specified direction
 *
 * @param type BC3D type
 * @param dir direction
 * @return bitmask
 */
#define BC3D_MASK(type, dir) ( (type) << (((dir)-1) * 3) )

#define BND3D_ID_ALL 0xFF0000000000000000 ///< BC3D: boundary id all
#define BOUND3D_ID(i) ((i)<<64) ///< BC3D: set boundary id for lift/drag @param i boundary ID

#endif
