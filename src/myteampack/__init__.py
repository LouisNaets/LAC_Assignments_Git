# pylint:disable=line-too-long
"""Main classes and functions.
"""
from pathlib import Path

from lacbox.htc import HTCFile
import numpy as np


class MyHTC(HTCFile):
    """Your team's class to auto-generate htc files.

    Instantiate a MyHTC object by passing in path to an htc file. E.g.:
        >> htc = MyHTC('./test.htc')
    """

    def _add_ctrltune_block(self, partial_load=(0.05, 0.7),  # pylint:disable=unused-argument
                            full_load=(0.06, 0.7),
                            gain_scheduling=2, constant_power=1,
                            regions=(1, 10, 12, 31), **kwargs):
        """Add a controller-tuning block to HAWC2S file. See HAWCStab2
        manual for explanation of input parameters.
        """

        ctrltune = self.hawcstab2.add_section('controller_tuning',
                                              pre_comments='\n; Need this block to call compute_controller_input')
        ctrltune.add_line(name='partial_load', values=partial_load, comments='fn [hz], zeta [-]')
        ctrltune.add_line(name='full_load', values=full_load, comments='fn [hz], zeta [-]')
        ctrltune.add_line(name='gain_scheduling', values=[gain_scheduling], comments='1 linear, 2 quadratic')
        ctrltune.add_line(name='constant_power', values=[constant_power], comments='0 constant torque, 1 constant power at full load')
        ctrltune.add_line(name='regions', values=regions, comments='Index of opt point (starting from 1) where new ctrl region starts')

    def _add_hawc2s_commands(self, rigid, gradient, tipcorr=True, induction=True,
                             num_modes=12, **kwargs):
        """Add commands for HAWC2S to end of hawcstab2 block.
        
        Comment/uncomment a HAWC2S command by passing the command as a
        keyword argument equal to True. E.g.:
            >> htc._add_hawc2s_commands(rigid=True, compute_steady_states=True, save_power=True)
        """
        defl_key = ['', 'no'][rigid]  # blade deformation or no?
        tipcorr_key = ['no', ''][tipcorr]  # tip correction or no?
        ind_key = ['no', ''][induction]  # induction or no?
        grad_key = ['no', ''][gradient]  # gradients or no?
        # add the pre-amble and output folder
        self.hawcstab2.add_line(name='', values=[''], comments='\n; HAWC2S commands (uncomment as needed)')
        self.hawcstab2.add_line(name='output_folder', values=['res_hawc2s'], comments='define the folder where generated files should be saved')
        # commands, options, and comments for hawc2s commands
        hawc2s_commands = {'compute_optimal_pitch_angle': [['use_operational_data'],
                                                           're-calculate and save opt file (pitch/rotor speed curve)'],
                           'compute_steady_states': [[f'{defl_key}bladedeform', f'{tipcorr_key}tipcorrect',
                                                      f'{ind_key}induction', f'{grad_key}gradients'],
                                                     'compute steady states -- needed for aeroelastic calculations'],
                           'save_power': [[''], 'save steady-state values to .pwr'],
                           'save_induction': [[''], 'save steady-state spanwise values to .ind files, 3 for each wind speed'],
                           'compute_stability_analysis': [[f'windturbine {num_modes}'],
                                                          'compute/save aeroelastic campbell diagram (.cmb), XX modes'],
                           'save_modal_amplitude': [[''], 'save modal amplitudes and phrases to .amb file'],
                           'compute_controller_input': [[''], 'calculate/save controller parameters (reqs. steady_states)'],
                           }
        # iterate over possible hawc2s commands
        for command, (values_, comments_) in hawc2s_commands.items():
            # if the command was passed in as a keyword argument, use the value passed in
            if command in kwargs:
                uncomment_command = kwargs[command]
            # if not given as keyword argument, default is for command to be commented
            else:
                uncomment_command = False
            name_ = [';', ''][uncomment_command] + command
            self.hawcstab2.add_line(name=name_, values=values_, comments=comments_)

    def _check_hawcstab2(self):
        """Verify HAWCStab2 block exists."""
        try:  # try to access the block
            self.hawcstab2
        except KeyError as exc:  # if we get a KeyError, block missing
            print('HAWCStab2 block not present in base file! Halting.')
            raise exc

    def _del_not_h2s_blocks(self):
        """Delete blocks in the htc file that HAWC2S doesn't use."""
        del self.simulation
        del self.dll
        del self.output
        del self.wind  # remove entire wind block for simplicity
        del self.aerodrag  # hawcstab2 dies with aerodrag can't handle this either

    def _set_initial_rotor_speed(self, omega0, bodyname='shaft', rotvec=(0., 0., -1.)):
        """Set the initial rotor speed on the shaft [rad/s]"""
        body = self.new_htc_structure.orientation.get_section(bodyname, field='mbdy2')
        body.mbdy2_ini_rotvec_d1 = rotvec + (omega0,)

    def _update_name_and_save(self, save_dir, append):
        """Update filename and save the file."""
        save_dir = Path(save_dir)  # sanitize inputs
        # get new name (excl extension) from HTCFile attribute "filename"
        name = Path(self.filename).name.replace('.htc', append)
        # set filename using HTCFile method
        self.set_name(name)
        # save the file
        self.save((save_dir / (name + '.htc')).as_posix())

    def make_hawc2s(self, save_dir, rigid, append, opt_path,
                    genspeed=(0, 480), minpitch=(0), opt_lambda=(8), **kwargs):
        """Make a HAWC2S file with specific settings.

        Args:
            save_dir (str/pathlib.Path): Path to folder where the htc file
                should be saved.
            rigid (boolean): Whether HAWC2S analysis should be a rigid or flexible
                structure.
            append (str): Text to append to the name of the master file.
            opt_path (str): Relative path from the saved htc file to the opt_file.
            genspeed (tuple, optional): 2-element tuple of minimum and maximum generator
                speed. Defaults to (0, 480).
        """
        # verify the file has hawcstab2 block
        self._check_hawcstab2()
        # delete blocks in master htc file that HAWC2S doesn't use
        self._del_not_h2s_blocks()
        # update the flexibility parameter in operational_data subblock
        defl_flag = [1, 0][rigid]  # 0 if rigid=True, else 1
        self.hawcstab2.operational_data.include_torsiondeform = defl_flag
        # correct the path to the opt file
        self.hawcstab2.operational_data_filename = opt_path
        # update the minimum generator speed
        self.hawcstab2.operational_data.genspeed = genspeed
        # update the minimum pitch
        self.hawcstab2.operational_data.minpitch = minpitch
        # update the tsr
        self.hawcstab2.operational_data.opt_lambda = opt_lambda
        # add hawc2s commands
        self._add_hawc2s_commands(rigid=rigid, **kwargs)
        # update filename and save the file
        self._update_name_and_save(save_dir, append)
        print(f'File "{append}" saved.')

    def make_hawc2s_ctrltune(self, save_dir, rigid, gradient, append, opt_path,
                    genspeed=(0, 480), minpitch=(0), opt_lambda=(8), gearratio=(1),
                    partial_load=(0.05, 0.7), full_load=(0.06, 0.7),
                    gain_scheduling=2, constant_power=1, regions=(1, 10, 12, 31), **kwargs):
        """Make a HAWC2S file with specific settings.

        Args:
            save_dir (str/pathlib.Path): Path to folder where the htc file
                should be saved.
            rigid (boolean): Whether HAWC2S analysis should be a rigid or flexible
                structure.
            append (str): Text to append to the name of the master file.
            opt_path (str): Relative path from the saved htc file to the opt_file.
            genspeed (tuple, optional): 2-element tuple of minimum and maximum generator
                speed. Defaults to (0, 480).
            partial_load (tuple, optional): Partial load parameters. Defaults to (0.05, 0.7).
            full_load (tuple, optional): Full load parameters. Defaults to (0.06, 0.7).
            gain_scheduling (int, optional): Gain scheduling parameter. Defaults to 2.
            constant_power (int, optional): Constant power parameter. Defaults to 1.
            regions (tuple, optional): Regions parameter. Defaults to (1, 10, 12, 31).
        """
        # verify the file has hawcstab2 block
        self._check_hawcstab2()
        # delete blocks in master htc file that HAWC2S doesn't use
        self._del_not_h2s_blocks()
        # update the flexibility parameter in operational_data subblock
        defl_flag = [1, 0][rigid]  # 0 if rigid=True, else 1
        self.hawcstab2.operational_data.include_torsiondeform = defl_flag
        # correct the path to the opt file
        self.hawcstab2.operational_data_filename = opt_path
        # update the minimum generator speed
        self.hawcstab2.operational_data.genspeed = genspeed
        # update the minimum pitch
        self.hawcstab2.operational_data.minpitch = minpitch
        # update the gear ratio
        self.hawcstab2.operational_data.gearratio = gearratio
        # update the tsr
        self.hawcstab2.operational_data.opt_lambda = opt_lambda
        # add hawc2s commands
        self._add_hawc2s_commands(rigid=rigid, gradient=gradient, **kwargs)
        # add controller tuning block
        self._add_ctrltune_block(partial_load=partial_load, full_load=full_load,
                                 gain_scheduling=gain_scheduling, constant_power=constant_power,
                                 regions=regions, **kwargs)
        # update filename and save the file
        self._update_name_and_save(save_dir, append)
        print(f'File "{append}" saved.')

    def _update_ctrl_params(self, cp_dict:dict):
        cp = self.dll.type2_dll__1.init # control param 
        cp.constant__2 = [2, 0] # min rotor compute
        cp.constant__3 = [3, 0.903] # rated rotor 8.627 RPM
        cp.constant__4 = [4, 10000/0.903*1.5*1000] # max gen torque multiplied by 1.5 for safety margin
        cp.constant__5 = [5, 0] # min pitch
        cp.constant__11 = [11, cp_dict["K_Nm/(rad/s)^2"]]
        cp.constant__12 = [12, cp_dict["KpTrq_Nm/(rad/s)"]]
        cp.constant__13 = [13, cp_dict["KiTrq_Nm/rad"]]
        cp.constant__16 = [16, cp_dict["KpPit_rad/(rad/s)"]]
        cp.constant__17 = [17, cp_dict["KiPit_rad/rad"]]
        cp.constant__21 = [21, cp_dict["K1_deg"]]
        cp.constant__22 = [22, cp_dict["K2_deg^2"]]

    def make_step(self, save_dir,  append, cp_dict, t_start:float, t_end:float, start_wsp:float, tint:float, turb_format:int, 
                  shear_format, tower_shadow_method:int, wind_ramp_abs):
        
        #update control parameterse
        self._update_ctrl_params(cp_dict)
        del self.hawcstab2
        self.simulation.time_stop = (t_end-t_start)+100
        self.wind.wsp = start_wsp
        self.wind.tint = tint
        self.wind.shear_format = shear_format
        self.wind.turb_format = turb_format 
        self.wind.tower_shadow_method = tower_shadow_method
        count = 0
        for v_o in range(wind_ramp_abs[2], wind_ramp_abs[3]):
            i = v_o-wind_ramp_abs[2]+1
            self.wind.add_line('wind_ramp_abs', values=(99+41*i, 100+41*i, 0, 1), comments=f'wsp after step {v_o+1}')
            count = count+1
        last_i = (wind_ramp_abs[3] - wind_ramp_abs[2] + 1)
        for v_o in range(wind_ramp_abs[3], wind_ramp_abs[2], -1):  # Go backwards
            i = last_i + (wind_ramp_abs[3] - v_o)  # Continue counting up from the ramp-up
            self.wind.add_line('wind_ramp_abs', values=(99 + 41 * i, 100 + 41 * i, 0, -1), comments=f'wsp after step {v_o - 1}')

        self._update_name_and_save(save_dir, append)
        print(f'File "{append}" saved.')