#! python3
import rawsHandler as hdl


class SensorClass(object):
    def __init__(self, channel, status, low_limit_of_measurement, high_limit_of_measurement, position):
        self.channel = channel
        self.status = status
        self.llm = low_limit_of_measurement
        self.hlm = high_limit_of_measurement
        self.position = position


if __name__ == "__main__":
    sensor_info = (hdl.opener(r'./sensor_status')).T
    sensors = list()
    for i in range(len(sensor_info)):
        sensors.append(SensorClass(*sensor_info[i]))


